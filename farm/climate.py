#%% Climate Class Definition
from math import exp
import numpy as np
from numpy.random import exponential, uniform

default_climate = {
    'alpha_r': [10.0] * 12,
    'lambda_r': [0.25] * 12,
    'doy_start': 1,
    't_sim': 180,
    'ET_max': 6.5
}

from datetime import timedelta, datetime
datetimes = np.arange(
    datetime(2018,1,1), datetime(2019,1,1), timedelta(days=1)
    ).astype(datetime)
month_value_by_day = np.array([datetime.month for datetime in datetimes])

class Climate():

    """ Creates a daily rainfall timeseries for use in ecohydrological modeling

    Usage: climate = Climate(alpha_r, lambda_r, t_seas, ET_max)

        alpha_r = average storm depth [mm]
        lambda_r = storm frequency [day^-1]
        t_sim = length of simulation [days]

    Default values:
        alpha_r = 10
        lambda_r = 0.25
        t_seas = 180
        ET_max = 6.5

    Note: lambda must either be a single value (constant rainfall probability all season),
    or have length of tseas (discrete rainfall probabilities each day.

    """
    def __init__(self, climate_parameters=default_climate, **kwargs):
        
       
        # Unpack the dictionary:
        self.t_sim = climate_parameters['t_sim'] 
        self.doy_start = climate_parameters['doy_start']
        self.ET_max = climate_parameters['ET_max']
       
        # Force doy to be in [1,365]:
        doy = np.arange(self.doy_start, self.doy_start+self.t_sim)
        while (doy - 365 > 0).any() == True:
            doy = doy - 365 * ((doy - 365) > 0)
        self.doy = doy

        # Determine daily values of lambda and alpha:
        lambda_r = climate_parameters['lambda_r']
        alpha_r = climate_parameters['alpha_r']
        # Check to ensure that lambda_r is either:
        # 1. a scalar, which means we have a constant climate
        # 2. has length of t_sim, which means we have specified daily values.
        # 3. has length of 12, which means we have specified monthly values.
        if isinstance(lambda_r, (float, int)):
            if isinstance(alpha_r, (float, int)):
                # We have a constant value:
                lambda_r_list = [lambda_r] * self.t_sim
                alpha_r_list = [alpha_r] * self.t_sim
            else:
                raise ValueError("lambda_r values and alpha_r values must be same length")
        elif len(lambda_r) == self.t_sim:
            if len(alpha_r) == self.t_sim:
                # We have daily values:
                lambda_r_list = lambda_r
                alpha_r_list = alpha_r
            else:
                raise ValueError("lambda_r values and alpha_r values must be same length")
        elif len(lambda_r) == 12:
            if len(alpha_r) == 12:
                # We have monthly values (remember that python is zero-indexed)
                lambda_all_year = np.array(
                    [lambda_r[month_value-1] for month_value in month_value_by_day]
                    )
                alpha_all_year = np.array(
                    [alpha_r[month_value-1] for month_value in month_value_by_day]
                    )
                lambda_r_list = lambda_all_year[self.doy-1]
                alpha_r_list = alpha_all_year[self.doy-1]
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
        self.rainfall = self.generate(self.alpha_r, self.lambda_r, self.t_sim, self.doy_start)

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
    def generate(alpha_r, lambda_r, t_sim, doy_start):
        """ Makes a time series of rainfall based on parameters

        Usage:
            generate(alpha_r, lambda_r, t_seas)

            alpha_r = average storm depth [mm]
            lambda_r = storm frequency [day^-1]
            t_seas = lenght of growing season [days]

        Note: lambda must either be a single value (constant rainfall probability all season),
        or have length of tseas (discrete rainfall probabilities each day.

        """

        amounts = [exponential(scale=alpha_value, size=1)[0] for alpha_value in alpha_r]
        rain_days = (uniform(low=0, high=1, size=t_sim) <= lambda_r).astype(int)
        return amounts * rain_days

