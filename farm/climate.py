#%% Climate Class Definition
from math import exp
import numpy as np
import pandas as pd
from numpy.random import exponential, uniform
from dateutil.relativedelta import *
import scipy.stats as st

default_climate = {
    'alpha_r': [10.0] * 12,
    'lambda_r': [0.25] * 12,
    'ET_max': 6.5
}

from datetime import timedelta, datetime

datetimes = np.arange(
    datetime(2018,1,1), datetime(2019,1,1), timedelta(days=1)
    ).astype(datetime)

month_value_by_day = np.array([datetime.month for datetime in datetimes])
week_value_by_day = np.array([datetime.isocalendar()[1] for datetime in datetimes])
dekad_value_by_day = np.array([datetime.timetuple().tm_yday//10+1 for datetime in datetimes])
# In the previous line of code, the first dekad is not repeated 10  times, so this works as a workaround:
dekad_value_by_day=np.insert(dekad_value_by_day,0,1)
dekad_value_by_day=np.delete(dekad_value_by_day,365)
# TODO Add semi_month_value_by_day if using 


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


def make_climate_parameters(
        station='OL JOGI FARM',
        data_file="data/CETRAD/CETRAD_rainfall.csv",
        year_min=30,
        interval='dekad'):

    """ Defines function that takes a rainfall station time series and returns alpha and lambda values by 
    a certain interval between week (7-days), dekad (10-days), semi-monthly (twice per month) or monthly.

        Usage:

            make_climate_parameters(
                station='OL JOGI FARM', 
                data_file="data/CETRAD/CETRAD_rainfall.csv",
                year_min= 30,
                interval='dekad' 
            )

        Default values:
            station = 'OL JOGI FARM' [string] # Rainfall Climatology for Laikipia 
            data_file = "data/CETRAD/CETRAD_rainfall.csv" # Path to file
            year_min = 30 # Minimum number of years required in timeseries
            interval = 'dekad' # Interval to calculate alphas andlambdas

                returns alpha_values, lambda_values
    """
    # Prepare the CETRAD dataset.
    year_min = year_min # minimum number of years to consider for a valid climate record.

    df = pd.read_csv(data_file)  # Read in the raw csv data.

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
    df['Week'] = [dt.week for dt in df['Datetime']]
    df['Semi_Month'] = (df['Datetime'].dt.day
                          .gt((df['Datetime']+pd.tseries.offsets.MonthEnd()).dt.day//2) 
                          + df['Month']*2 -1)
    df['Dekad'] = df['Datetime'].dt.dayofyear//10+1
    
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

    columns = [station] + ['Year', 'Month', 'Week', 'Dekad', 'Semi_Month','Datetime']
    rainfall = df[columns]

    # First, find all the rows in the data where it rained and group by month.
    rain_days = rainfall.loc[rainfall[station] > 0]

    # Find all locations in the data where an observation was made.
    all_days = rainfall.loc[rainfall[station] >= 0] 

    # Find just the rainfall amounts on days that it rained.
    data = rainfall.loc[rainfall[station] > 0][station]
    
    # Fit the daily rainfall amounts to an exponential distribution.
    check_exponential(data)

    if interval == 'month':
        # Determine the Monthly values of alpha and lambda from the station data:
        lambda_values = (
            rain_days.groupby('Month')[station].count() /
            all_days.groupby('Month')[station].count()
        )
        alpha_values = rain_days.groupby('Month')[station].mean()
    elif interval == 'dekad':
        lambda_values = (
            rain_days.groupby('Dekad')[station].count() / 
            all_days.groupby('Dekad')[station].count()
        )
        alpha_values = rain_days.groupby('Dekad')[station].mean()
    elif interval == 'semi_month':
        lambda_values = (
            rain_days.groupby('Semi_Month')[station].count() / 
            all_days.groupby('Semi_Month')[station].count()
        )
    else:
        raise(NotImplementedError)

    return alpha_values.to_list(), lambda_values.to_list()


# Should we rename this?
class Climate():

    """ Creates a years worth of daily rainfall timeseries for use in ecohydrological modeling

    Usage: climate = Climate(alpha_r, lambda_r, ET_max)

        alpha_r = average storm depth [mm]
        lambda_r = storm frequency [day^-1]
      
    Default values:
        alpha_r = 10
        lambda_r = 0.25
        t_seas = 180
        ET_max = 6.5

    Note: lambda must either be a single value (constant rainfall probability all season),
    or have length of tseas (discrete rainfall probabilities each day.

    """
    def __init__(self, alpha_r=[10.0] * 12, lambda_r=[0.25] * 12, ET_max=6.5, **kwargs):
        
       

        self.ET_max = ET_max
        
        # First check to see if we were passed a station. If
        # so, then try to get the climate parameters from the station 
        # using the make_climate_parameters helper function.
        # NOTE: THIS WILL ONLY WORK FOR LAIKIPIA STATION DATA!
        #
        # Otherwise, check to ensure that lambda_r is either:
        # 1. a scalar, which means we have constant climate parameters
        # 2. has length of 365, which means we have specified daily values.
        # 3. has length of 12, which means we have specified monthly values.
        if 'station' in kwargs:
            if 'data_file' in kwargs:
                alpha_r, lambda_r = make_climate_parameters(
                    station=kwargs['station'],
                    data_file=kwargs['data_file'],
                    interval=kwargs['interval'])
            else:
                alpha_r, lambda_r = make_climate_parameters(station=kwargs['station'], interval=kwargs['interval'])
        if isinstance(lambda_r, (float, int)):
            if isinstance(alpha_r, (float, int)):
                # We have a constant value:
                lambda_r_list = [lambda_r] * 365
                alpha_r_list = [alpha_r] * 365
            else:
                raise ValueError("lambda_r values and alpha_r values must be same length")
        elif len(lambda_r) == 365:
            if len(alpha_r) == 365:
                lambda_r_list = lambda_r
                alpha_r_list = alpha_r
            else:
                raise ValueError("lambda_r values and alpha_r values must be same length")
        elif len(lambda_r) == 12:
            if len(alpha_r) == 12:
                # We have monthly values (remember that python is zero-indexed)
                lambda_r_list = np.array(
                    [lambda_r[month_value-1] for month_value in month_value_by_day]
                    )
                alpha_r_list = np.array(
                    [alpha_r[month_value-1] for month_value in month_value_by_day]
                    )
            else:
                raise ValueError("lambda_r values and alpha_r values must be same length")
        elif len(lambda_r) == 53:
            if len(alpha_r) == 53:
                # We have weekly values (remember that python is zero-indexed)
                lambda_r_list = np.array(
                    [lambda_r[week_value-1] for week_value in week_value_by_day]
                    )
                alpha_r_list = np.array(
                    [alpha_r[week_value-1] for week_value in week_value_by_day]
                    )
            else:
                raise ValueError("lambda_r values and alpha_r values must be same length")
        elif len(lambda_r) == 37:
            if len(alpha_r) == 37:
                # We have dekadal values (remember that python is zero-indexed)
                lambda_r_list = np.array(
                    [lambda_r[dekad_value-1] for dekad_value in dekad_value_by_day]
                    )
                alpha_r_list = np.array(
                    [alpha_r[dekad_value-1] for dekad_value in dekad_value_by_day]
                    )
            else:
                raise ValueError("lambda_r values and alpha_r values must be same length")
        else:
            raise ValueError(
                "lambda_r & alpha_r values should be a constant, have length of t_sim, or have length of 12"
            )

        # Set the rainfall parameters for this instance
        self.alpha_r = alpha_r_list
        self.lambda_r = lambda_r_list
        
        # Use the static method, generate, to create this instance's rainfall.
        self.rainfall = self.generate(self.alpha_r, self.lambda_r)

        # Assign any other passed parameters (e.g. site, etc...)
        for key, value in kwargs.items():
            self.key = value


    def calc_E(self, s, q=1.5, LAI=None, sh=None): 
        """ Determines the daily evaporation as a function of relative soil moisture

        Usage: calc_E(s)

            TODO: Update parameters.

            s = relative soil moisture [0-1]
            E_max_p = E_max * exp(-k * LAI)
            E = E_max * [(s-sh)/(1-sh)]^q

        """
        if LAI == None:
            raise ValueError("Climate calc_E expects LAI that's not None.")

        k = 0.5
        E_max_p = self.ET_max*exp(-k*LAI) 
        if s >= sh:
            return pow((s-sh)/(1-sh), q)*E_max_p
        else:
            return 0

    @staticmethod # Static methods can be called without instancing the class.
    def generate(alpha_r, lambda_r, t_sim=365, doy_start=1):
        """ Makes a time series of rainfall based on parameters

        Usage:
            generate(alpha_r, lambda_r, t_seas)

            alpha_r = average storm depth [mm]
            lambda_r = storm frequency [day^-1]
            t_seas = lenght of growing season [days]

        Note: lambda must either be a single value (constant rainfall probability all season),
        or have length of tseas (discrete rainfall probabilities each day.

        """
        # Force doy to be in [1,365]:
        doys = np.arange(doy_start, doy_start + t_sim)
        while (doys - 365 > 0).any() == True:
            doys = doys - 365 * ((doys - 365) > 0)

        amounts = [exponential(scale=alpha_r[doy-1], size=1)[0] for doy in doys]
        rain_days = [(uniform(low=0, high=1, size=1) <= lambda_r[doy-1] ).astype(int) for doy in doys]
        return np.multiply(amounts, [v[0] for v in rain_days])

