# Creates climate parameters for the model based on CETRAD rainfall data.

# List of possible stations:

"""
['ARCHERS POST', 'ARDENCAPLE FARM', 'CASTLE FOREST STN', 'CHOGORIA FOREST STN', 'CHUKA FOREST STN', 
'COLCHECCIO', 'DOL DOL DAO', 'EL KARAMA', 'EMBORI FARM', 'EMBU MET STN', 'ENASOIT FARM', 
'GATHIURU FOREST STN', 'HOMBE FOREST STN', 'IRANGI FOREST STN', 'ISIOLO DAO', 'JACOBSON FARM', 
'JUNCTION (EWASO NAROK)', 'KABARU FOREST STN', 'KAGURU', 'KALALU (NRM)', 'KAMWAKI FARM', 
'KARURI (NRM)', 'KINAMBA MOW', 'KISIMA FARM', 'LAMURIA MET STN', 'LARIAK FOREST STN', 
'LOGILADO (NRM)', 'LOLDAIGA FARM', 'LOLDOTO FARM', 'LOLMARIK FARM', 'LORUKU FARM', 
'MARALAL DC', 'MARIENE CRS', 'MATANYA (NRM)', 'MERU FOREST STN', 'MOGWONI RANCH', 'MPALA FARM',
 'MUGIE RANCH', 'MUKENYA FARM', 'MUKOGODO (NRM)', 'MUNYAKA (NRM)', 'MURINGATO FOREST STN', 
 'MUTARA ADC FARM', 'MWEA IRRIGATION SCHEME', 'NANYUKI FOREST STN', 'NANYUKI KAF', 
 'NARO MORU FG POST', 'NARO MORU FOREST STN', 'NARO MORU GATE STN', 'NARO MORU MET STN', 
 'NDARAGWA FOREST STN', 'NGENIA (NRM)', 'NGENIA B', 'NICOLSON FARM', 'NYERI MOW', 
 'OL ARABEL FOREST STN', 'OL BOLOSAT FOREST STN', 'OL DONYO FARM', 'OL JOGI FARM', 'OL JORO OROK FTC',
  'OL MYSOR FARM', 'OL PEJETA FARM', 'ONTULILI FOREST STN', 'PYRAMID OL JOGI', 
  'RAGATI FOREST STN', 'RUMURUTI (NRM)', 'RUMURUTI MOW', 'SATIMA FARM', 'SEGERA PLANTATIONS', 
  'SHAMATA', 'SIRAJI (NRM)', 'SIRIMA (NRM)', 'SOLIO RANCH', 'SOUTH MARMANET FOREST STN', 
  'SUGUROI ESTATE', 'TELEKI (MT KENYA)', 'TELESWANI (NRM)', 'THARUA FARM', 'TIMAU MARANIA', 'TRENCH FARM']
"""
import scipy.stats as st
import pandas as pd
import numpy as np
import copy
from datetime import datetime
from dateutil.relativedelta import *
import functools
from .climate import Climate
from .model import CropModel

def check_exponential(data):

	""" Defines function that fits daily rainfall amounts to an exponential distribution and returns pdf 
		and r2. The r2 should be above 0.9 to be an exponential.

		Usage:

			check_exponential(data):

				returns r2, pdf

		How it works:
		- Step 1: To fit the distribution, we use functions from python's suite of numerical analysis, scipy.
		The scipy.stats module has a large suite of distribution functions pre-defined, which we can use to 
		develop a fit for our data. The distribution we are interested in is the exponential distribution, 
		which is called expon in the stats module.

		- Step 2-4: Calculate fitted PDF and error with fit in distribution. To test the fit of our distribution, 
		we can compare the empirical histogram to that predicted by our model. We first use our `data` to generate 
		the empirical histogram. In this example, we break the data into `30` bins, and we generate a histrogram 
		of `density` rather than counts. This allows for an easier comparison between our empirical data and the 
		fitted probability distribution function. 
		
		Here are the steps:

		1. Generate a histogram, from the `data`. Save the bin locations in `x` and the density of values in `y`
		2. Shift the `x` bin locations generated from the histogram to the center of bins.
		3. Calculate the value of the fitted `pdf(x)` for each of the bins in `x`.
		4. Determine the residual sum of the squares, $SS_{error}$, and total sum of squares, $SS_{yy}$, according 
		to the equations in rainfall-variability.ipynb.
	"""

	# Step 1. Fit the distribution.
	distribution = st.expon
	params = distribution.fit(data, loc=0) # Force the distribution to be built off of zero

	arg = params[:-2]
	loc = params[-2]
	scale = params[-1]

	y, x = np.histogram(data, bins=30, density=True)

	# Step 2. Shift the x bin locations to the center of bins.
	x = (x + np.roll(x, -1))[:-1] / 2.0

	# Step 3. Calculate the values of pdx(x) for all x.
	pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)

	# Step 4. Determine the residual and total sum of the squares.
	ss_error = np.sum(np.power(y - pdf, 2.0))
	ss_yy = np.sum(np.power(y - y.mean(), 2.0))

	r_2 = 1 - ( ss_error / ss_yy )

	if r_2 < 0.9:
		print("WARNING. r2 for {station} is {r_2}".format(
			station=station,
			r_2=r_2))

	return r_2, pdf


def make_climate_parameters(station='OL JOGI FARM'):

	# Prepare the CETRAD dataset.
	year_min = 30 # minimum number of years to consider for a valid climate record.

	df = pd.read_csv("../data/CETRAD/CETRAD_rainfall.csv")  # Read in the raw csv data.

	# Step 1. Convert text strings into datetime objects.
	format = '%m/%d/%y' # Column RDate has data in M/D/YY
	df['Datetime']=pd.to_datetime(df['RDate'], format=format) # Create a new column of datetime objects using RDate.

	# 2. Step 2. Convert future dates inferred during the conversion back into 20th century dates.
	# Python is a future-looking programming language, and assumes that 1/1/34 is Jan 1, 2034.
	# We can fix this by finding all the dates in the future (dt > datetime.now()) and removing 100 years from
	# their value. This requires using the relativedelta function, which handles weird stuff like leap years.
	df['Datetime'] = df['Datetime'].map(lambda dt: dt+relativedelta(years=-100) if dt > datetime.now() else dt)

	# Step 3. Extract the Year and Month from the Datetime to make aggregation easier.
	df['Year'] = [dt.year for dt in df['Datetime']]
	df['Month'] = [dt.month for dt in df['Datetime']]

	n_years = len(df['Year'].unique())

	# Check to make sure we have enough data for fitting and parameter estimation.
	if n_years < year_min:
		print("WARNING! Station record for {station} has only {n_years} years.".format(
			station=station,
			n_years=n_years))

	# Step 4. Use the Datetime values as the index for this dataframe.
	df = df.set_index(pd.DatetimeIndex(df['Datetime']))  # Set the Datetime column as the dataframe index

	# Step 5.  Delete the old RDate column, which we no longer need. 
	# We will keep the Datetime column, in case we need it later.
	df = df.drop(['RDate'], axis=1)

	columns = [station] + ['Year', 'Month', 'Datetime']
	rainfall = df[columns]

	# First, find all the rows in the data where it rained and group by month.
	rain_days = rainfall.loc[rainfall[station] > 0]

	# Find all locations in the data where an observation was made.
	all_days = rainfall.loc[rainfall[station] >= 0]

	# Find just the rainfall amounts on days that it rained.
	data = rainfall.loc[rainfall[station] > 0][station]
	
	# Fit the daily rainfall amounts to an exponential distribution.
	check_exponential(data)

	# Determine the Monthly values of alpha and lambda from the station data:
	lambda_by_month = (
	    rain_days.groupby('Month')[station].count() /
	    all_days.groupby('Month')[station].count()
	)

	alpha_by_month = rain_days.groupby('Month')[station].mean()

	# MAKE THE CLIMATE PARAMETER DICT:
	climate = pd.DataFrame(alpha_by_month)
	climate = climate.rename(columns={'OL JOGI FARM': 'alpha_by_month'})
	climate['lambda_by_month'] = lambda_by_month
	

	return climate['alpha_by_month'].to_list(), climate['lambda_by_month'].to_list()


@functools.lru_cache(maxsize=128)
def average_soil_moisture(model, n_sims=100, t_before=30, doy=None):

	alpha_r = model.climate.alpha_r
	lambda_r = model.climate.lambda_r
	climates = [Climate(alpha_r, lambda_r) for sim in np.arange(n_sims)]
    
	# Create a temporary crop object with a 0 day length of growing period.
	temp_crop = copy.copy(model.crop)
	temp_crop.lgp = 0

    # Get output from each simulataion using an implicit for loop.
	# Use the temp crop object to create these models.
	models = [ CropModel(crop=temp_crop,soil=model.soil,climate=climates[i]) for i in np.arange(n_sims) ]
	
	output = [ models[i].run(do_output=True, planting_date=doy+1, t_before=t_before, t_after=0) for i in np.arange(n_sims) ]

    # Extract the final value of soil moisture from each output.
	values = pd.DataFrame([output[i]['s'][-1:] for i in np.arange(n_sims)])
	return values.mean(), values.std()