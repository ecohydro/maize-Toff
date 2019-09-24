#%% Climate Class Definition
from math import exp
from numpy.random import exponential, uniform

class Climate():
    """ Creates a daily rainfall timeseries for use in ecohydrological modeling

    Usage: climate = Climate(alpha_r,lambda_r,t_seas,ET_max)

        alpha_r = average storm depth [mm]
        lambda_r = storm frequency [day^-1]
        t_seas = length of growing season [days]

    Default values:
        alpha_r = 10
        lambda_r = 0.25
        t_seas = 180
        ET_max = 6.5

    Note: lambda must either be a single value (constant rainfall probability all season),
    or have length of tseas (discrete rainfall probabilities each day.

    """
    def __init__(self, alpha_r=10.0, lambda_r=0.25, t_seas=180, ET_max=6.5, **kwargs):
        # Check to ensure that lambda_r is either a scalar or has length of t_seas
        if not isinstance(lambda_r, float):
            if len(lambda_r) != t_seas:
                raise ValueError("lamda_r values should be a constant, or have length of t_seas")
        # Set the rainfall parameters for this instance
        self.alpha_r = alpha_r
        self.lambda_r = lambda_r
        self.t_seas = t_seas
        self.ET_max = ET_max
        # Use the static method, generate, to create this instance's rainfall.
        self.rainfall = self.generate(self.alpha_r, self.lambda_r, self.t_seas)

        # Assign any other passed parameters (e.g. site, etc...)
        for key, value in kwargs.items():
            self.key = value


    def calc_E(self, s, t, q=4, sh=None, LAI=None):  # soil=None, plant=None): 
        """ Determines the daily evaporation as a function of relative soil moisture

        Usage: calc_E(s)

            TODO: Update parameters.

            s = relative soil moisture [0-1]
            t = day of season

            E_max_p = E_max * exp(-k * LAI)
            E = E_max * [(s-sh)/(1-sh)]^q

        """
        if LAI == None:
            raise ValueError("Climate calc_E expects LAI that's not None.")

        k = -0.5
        #print("k:",-k)
        #print("ET_max:",self.ET_max)
        #print("LAI:",LAI)
        E_max_p = self.ET_max*exp(-k*LAI) # plant.calc_LAI(t)/plant.LAI_max)
        #print(E_max_p)

        return pow((s-sh)/(1-sh), q)*E_max_p

    @staticmethod # Static methods can be called without instancing the class.
    def generate(alpha_r, lambda_r, t_seas):
        """ Makes a time series of rainfall based on parameters

        Usage:
            generate(alpha_r, lambda_r, t_seas)

            alpha_r = average storm depth [mm]
            lambda_r = storm frequency [day^-1]
            t_seas = lenght of growing season [days]

        Note: lambda must either be a single value (constant rainfall probability all season),
        or have length of tseas (discrete rainfall probabilities each day.

        """
        amounts = exponential(scale=alpha_r, size=t_seas)
        rain_days = (uniform(low=0, high=1, size=t_seas) <= lambda_r).astype(int)
        return amounts * rain_days
