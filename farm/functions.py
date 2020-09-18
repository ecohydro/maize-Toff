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
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import copy
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import proplot as plot
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
	); print(lambda_by_month)

	alpha_by_month = rain_days.groupby('Month')[station].mean()

	# MAKE THE CLIMATE PARAMETER DICT:
	climate = pd.DataFrame(alpha_by_month)
	climate = climate.rename(columns={station: 'alpha_by_month'})
	climate['lambda_by_month'] = lambda_by_month

	return climate['alpha_by_month'].to_list(), climate['lambda_by_month'].to_list(), rainfall


@functools.lru_cache(maxsize=128)
def average_soil_moisture(model, n_sims=100, t_before=60, doy=None):

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

def calc_yield(stress=None, max_yield = 4680):
    yield_kg_ha = -max_yield*stress + max_yield
    
    if stress > 1:
        raise ValueError("static stress, {stress} is larger than 1".format(
                stress=stress))
    if stress < 0:
        raise ValueError("static stress, {stress} is less than 0".format(
                stress=stress))
    
    return yield_kg_ha

# TODO: Consider moving plotting functions into their own script.
def plot_lin_regression(x_var = None, y_var = None, x_str = None, y_str = None, data = None, 
                        ann_x = 101, ann_y = 4500, 
                        x_lab = 'X label here', y_lab = 'Y label here', title = 'Title here', positive = True, plot = False):
    """ Computes linear regression between independent and dependent variable. 
    Usage: plot_lin_regression(x_var, y_var, x_lab, y_lab, title)
        ann_x = where on x-axis annotation should be placed
        ann_y = where on y-axis annotation should be placed
        Returns: R_squared, m, b
    """
    # Define variables
    X, y = x_var, y_var
    
    # Linear regression
    denominator = X.dot(X) - X.mean() * X.sum()
    m = ( X.dot(y) - y.mean() * X.sum() ) / denominator
    b = (y.mean() * X.dot(X) - X.mean() * X.dot(y) ) / denominator

    y_pred = m*X + b

    if plot == False:
        # Calculate residuals
        res = y - y_pred
        tot = y - y.mean()

        R_squared = 1 - res.dot(res) / tot.dot(tot)

    else:
        plt.figure(figsize=(5,4))

        g = sns.lmplot(x_str, y_str, data, ci=95, height=4, scatter_kws={'color':'black','alpha':0.6}) # ,, line_kws={'color': 'black'}
     
        # Calculate residuals
        res = y - y_pred
        tot = y - y.mean()

        R_squared = 1 - res.dot(res) / tot.dot(tot)
        print(R_squared)
        print('m',m)
        print('b',b)

        props = dict(boxstyle='square', facecolor='white', alpha=0.5, lw = 1.5) # , ec="b"

        # place a text box in upper left in axes coords
        plt.text(ann_x, ann_y, textstr, fontsize=10, #transform=ax.transAxes, 
                verticalalignment='top', bbox=props)

        plt.xlabel(x_lab)
        plt.ylabel(y_lab)
        plt.title(title, fontweight="bold")
    
    if positive == True:
        textstr = '\n'.join((
            r'$ y = %.2f$x' % (m, )+'+$  %2.0f$' % (b, ),
            r'$r^2 = %.2f$' % (R_squared, ))) 
    else:
        textstr = '\n'.join((
        r'$ y = %.2f$x' % (m, )+'$  %2.0f$' % (b, ),
        r'$r^2 = %.2f$' % (R_squared, )))

    return R_squared, m, b

def power_law_fit(xdat,ydat, x_lab, y_lab, title):
    x,y = xdat, ydat
    #xmax = 4260
    xmax = max(x)
    power_law = lambda x, a, b: a * (x**b)
    
    f, axs = plot.subplots(ncols=1, nrows=2, share=0, figsize=(5,4)) 
    
    # Find best fit.
    popt, pcov = curve_fit(power_law, x, y)
    
    # Top plot
    # Plot data and best fit curve.
    axs[0].plot(x, y,'ok', alpha=0.6)
    axs[0].plot(np.sort(x), power_law(np.sort(x), *popt),'-',markersize=3,  linewidth=2.5) # like this color color=(0.2, 0.4, 0.6, 0.6)
    axs[0].format(xlim=(100, 1000), ylim=(0, 5000))

    #r2, v1
    residuals = y - power_law(x, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y-np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    r_squared
    
    #r2, v2
    from sklearn.metrics import r2_score
    r2_score(y, power_law(x, *popt), multioutput='variance_weighted')
    
    # Add text
    textstr = r'$r^2=%.2f$' % (r_squared, )
    props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
    axs[0].format(suptitle=title, title = textstr,titleweight='bold', titleloc='ul',
                 ylabel=y_lab, xlabel=x_lab)
    
    # Bottom plot
    axs[1].plot(residuals) #linewidth=.9
    axs[1].format(title='Residuals', titleweight='bold',xlabel='Simulation Number',
                 ylabel='Error') #, titleloc='ul
    axs[0].set_xlim(min(x)-3, max(x)+10)  


def plot_newfit(xdat,ydat, x_lab, y_lab, title=None):
    x,y = xdat, ydat
    #xmax = 4260
    x,y = xdat, ydat
    xmax = x.max()

    def new_fit(x, A, B):
        return A*(x - xmax)**2+B # testing this out

    f, axs = plot.subplots(ncols=1, nrows=3, share=0, figsize=(5,4)) 
    
    # Find best fit.
    popt, pcov = curve_fit(new_fit, x, y)
    
    # Top plot
    # Plot data and best fit curve.
    axs[0].plot(x, y,'ok', alpha=0.6)
    axs[0].plot(np.sort(x), new_fit(np.sort(x), *popt),'-',markersize=3,  linewidth=2.5) # like this color color=(0.2, 0.4, 0.6, 0.6)
    axs[0].format(xlim=(0, 1000), ylim=(0, 5000))
    
    #r2, v1
    residuals = y - new_fit(x, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y-np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    r_squared
    
    #r2, v2
    from sklearn.metrics import r2_score
    r2_score(y, new_fit(x, *popt), multioutput='variance_weighted')
    
    # Add text
    textstr = r'$r^2=%.2f$' % (r_squared, )
    props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
    axs[0].format(suptitle='y = A(x - xmax)**2B', title = textstr,titleweight='bold', titleloc='ul',
                 ylabel=y_lab, xlabel=x_lab)
    
    # Bottom plot
    axs[1].plot(residuals) #linewidth=.9
    axs[1].format(title='Residuals, Time Series', titleweight='bold',xlabel='Simulation Number',
                 ylabel='Error') #, titleloc='ul
    axs[0].set_xlim(min(x)-3, max(x)+10)  

    # Plot histogram of the residuals
    n_bins = 100
    axs[2].hist(residuals, bins=n_bins)
    axs[2].format(title='Residuals, Histogram', titleweight='bold',xlabel='Error',
                 ylabel='Number of Simulations') #, titleloc='ul
    
    return residuals
    

def polyfit(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)

     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot

    return results, yhat, ybar

def plot_polyfit(x=None, y=None, degree=None, x_lab='Seasonal rainfall (mm)',y_lab='Yield (kg/ha)',title='Polynomial fit'):
    # degree = degree of the fitting polynomial
    # stop the linspace at the x max of data
    xmin = min(x)
    xmax = max(x)

    fig, ax = plt.subplots(figsize=(5,4))
    p = np.poly1d(np.polyfit(x, y, degree))
    t = np.linspace(xmin, xmax, 1000)

    ax.plot(x, y, 'ok', t, p(t), '-', markersize=3, alpha=0.6, linewidth=2.5)
    #ax.format(xlim=(0, xmax+50), ylim=(0, 4000))
    results, yhat, ybar = polyfit(x,y,degree)
    

    R_squared = results['determination']
    textstr = r'$r^2=%.2f$' % (R_squared, )
    props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)

    fig.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)
    
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.title(title, fontweight="bold")

    results['polynomial'][0]

    # TODO: confidence intervals around line?


def evolved_calc_yield(dtm=None, m=None, b=None):
    yield_kg_ha = m*dtm + b

    if dtm > 185:
        raise ValueError("days to maturity, {dtm} is larger than 175".format(
                dtm=dtm))
    if dtm < 68:
        raise ValueError("days to maturity, {dtm} is less than 68".format(
                dtm=dtm))

    return yield_kg_ha

# verified using Kenya Seed Co. - https://web.archive.org/web/20190819125927/http://kenyaseed.com/gallery/maize/
verified_hybrid_data = pd.read_csv('../data/Yields/hybrid_yields_verified.csv')

# convert to metric tons
verified_hybrid_data['yield_metric_tons'] = verified_hybrid_data.verified_yield_kg_acre/1000
p, m, b = plot_lin_regression(verified_hybrid_data.verified_days_to_maturity, verified_hybrid_data.yield_metric_tons, 
                             'verified_days_to_maturity', 'yield_metric_tons', verified_hybrid_data, 
                             85, 4.9, 'Days to Maturity (days)', 'Yield (tons/ha)', 
                             'Potential Maize Yields from Kenya Seed Company', positive=False)