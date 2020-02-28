#%% Climate Class Definition
from math import exp
import numpy as np
from numpy.random import exponential, uniform

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
        
       
        # Unpack the dictionary:
        self.ET_max = ET_max
        
        # Check to ensure that lambda_r is either:
        # 1. a scalar, which means we have constant climate parameters
        # 2. has length of 365, which means we have specified daily values.
        # 3. has length of 12, which means we have specified monthly values.
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

