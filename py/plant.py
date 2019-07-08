#%% Defines a Plant class
from py.soil import Soil

#%%
class Plant():
    """ Defines a plant class


    """
    def __init__(self,
        Zr=500,             # Rooting depth [mm]
        T_max=4,            # Maximum transpiration [mm/day]
        sw_MPa = -1.5,      # wilting point of plant in water potential [MPa]
        s_star_MPa = -0.2,  # water potential of maximum transpiration [MPa]
        soil=None,          # a soil in which this plant will grow
        climate=None        # a climate in which this plant will grow 
    ):
        self.Zr = Zr
        self.sw_MPa = sw_MPa
        self.s_star_MPa = s_star_MPa

        # Use the soil to determine critical plant parameters
        self.sw = soil.s(soil.theta(sw_MPa))
        self.s_star = soil.s(soil.theta(s_star_MPa)) 

        self.T_max = T_max      # Should be over-written by any subclass.
        
    def calc_LAI(self):
        raise NotImplementedError
    
    def calc_T(self, s, climate):
        raise NotImplementedError

#%%
class Crop(Plant):
    """ Creates a Crop class.

    Usage: crop = Crop(
        kc_max=1,
        LAI_max=3.0,
        day_of_season=1,
        climate=climate,
        soil=soil
    )

    """
    def __init__(self,*args,**kwargs):
        self.kc_max = kwargs.pop('kc_max' )     # Maximum crop coefficient [0-1]
        self.LAI_max = kwargs.pop('LAI_max')    # Maximum crop leaf area index [m^2/m^2]
        climate = kwargs.pop('climate')         # The climate that this plant will grow in
        super(Crop, self).__init__(*args, **kwargs)
    
        self.kc = self.calc_kc(day_of_season)
        self.LAI = self.calc_LAI(self.kc)
        self.T_max = self.calc_T_max(climate)

    def calc_kc(self, day_of_season):
        """ Calculates the current crop coefficient based on day of season

        Usage: calc_kc(day_of_season)

        Note: Currently just returns kc_max regardless of day of season.

        """
        return self.kc_max
    
    def calc_T_max(self, climate):
        """ Calculates max Transpiration variable.
        
        Usage: calc_T_max(climate)
            
        """
        return self.kc * climate.ET_max

    def calc_LAI(self, p=1):
        """ Returns a Leaf Area Index (LAI) variable. LAI comes
            from function of kc. Currently based on linear relationship 
            between kc and LAI (assumption).
        
        Usage: calc_LAI(p=1)

        LAI = (LAI_max/kc_max)^p * kc

        Note: p=1 assumes a linear relationship between LAI and kc

        """
        return pow((self.LAI_max/self.kc_max),p) * self.kc

    def calc_T(self, s, climate):
        """ Calculates Transpiration variable as a stepwise
            linear function.
        
        Usage: calc_T(s, climate)
        
            s = relative soil moisture [0-1]
            climate = a climate object used to determine maximum transpiration
        
        """
        if not s <= 1:
            raise ValueError("s must be <= 1")
        if s>=self.s_star:
            return self.calc_T_max(climate)
        elif s>=self.sw:
            return (s-self.sw)/(self.s_star-self.sw)*self.calc_T_max(climate)
        else:
            return 0

