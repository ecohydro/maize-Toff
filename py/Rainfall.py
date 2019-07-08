#%%
# Climate class for ecohydrological model.
# 
#%% Climate Class Definition

class Rainfall():
    """ Creates a daily rainfall timeseries for use in ecohydrological modeling

    Usage: rainfall = Rainfall(alpha_r,lambda_r,t_seas)

        alpha_r = average storm depth [mm]
        lambda_r = storm frequency [day^-1]
        t_seas = length of growing season [days]
    
    Note: lambda must either be a single value (constant rainfall probability all season),
    or have length of tseas (discrete rainfall probabilities each day. 

    """
    def __init__(self,alpha_r=10.0, lambda_r=0.3, t_seas=180, **kwargs):
        # Check to ensure that lambda_r is either a scalar or has lenght of t_seas
        if not isinstance(lambda_r, float):
            if len(lambda_r) != t_seas:
                raise ValueError("lamda_r values should be a constant, or have length of t_seas")
        # Set the rainfall parameters for this instance
        self.alpha_r = alpha_r
        self.lambda_r = lambda_r
        self.t_seas = t_seas
        self.rainfall = self.generate(self.alpha_r, self.lambda_r, self.t_seas)

        # Assign any other passed parameters (e.g. site, etc...)
        for key, value in kwargs.items():
            self.key = value
    
        
    @staticmethod
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
        # Import necessary functions from numpy
        from numpy.random import exponential, uniform

        amounts = exponential(scale=self.alpha_r, size=t_seas)
        rain_days = (uniform(low=0, high=1, size=t_seas) <= lambda_r).astype(int)
        return amounts * rain_days


        