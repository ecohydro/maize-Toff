#%%
class Plant():
    """ Defines a plant class


    """
    def __init__(self,
        Zr=500,             # Rooting depth [mm]
        T_max=4,            # Maximum transpiration [mm/day]
        sw_MPa = -1.5,      # wilting point of plant in water potential [MPa]
        s_star_MPa = -0.2,  # water potential of maximum transpiration [MPa]
        soil=None           # a soil in which this plant will grow
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
    
    def calc_T(self, s):
        raise NotImplementedError

#%%
class Crop(Plant):
    """ Creates a Crop class.

    Usage: crop = Crop(soil=soil)

    optional keyword arguments and their default values:
        'Zr': 500,          # Planting depth [mm]
        'sw_MPa':-1.5,      # Plant wilting point [MPa]
        's_star_MPa':-0.2,  # Water potential for max T
        'kc_max':1.2,       # Maximum crop coefficient
        'LAI_max':3.0,      # Max Leaf Area Index [m2/m2]
        'T_max':4.0         # Max Crop Water Use [mm/day]

    """
    def __init__(self, Zr=500, sw_MPa=-1.5, s_star_MPa=-0.2, kc_max=1.2, LAI_max=3.0, T_max=4.0,*args,**kwargs):
        self.kc_max = kc_max     # Maximum crop coefficient [0-1]
        self.LAI_max = LAI_max    # Maximum crop leaf area index [m^2/m^2]
        self.T_max = T_max        # Maximum crop water use [mm/day]
        super(Crop, self).__init__(*args, **kwargs)

    def calc_kc(self, day_of_season, t_seas = 120, f1 = 0.2, f2 = 0.5, f3 = 0.75, EoS = 1.0, kc_ini = 0.30, kc_max = 1.2, kc_EoS = 0.6):
        """ Calculates crop coefficient that varies throughout the season 
        
        Usage: calc_kc(self, day_of_season, t_seas = 120, f1 = 0.2, f2 = 0.5, f3 = 0.75, EoS = 1.0, kc_ini = 0.30, kc_max = 1.2, kc_EoS = 0.6)
            Note: t must be a single-dimension array
            day_of_season = user input # Start date [day]
            t_seas = 120    # Length of growing season [days]
            f1 = 0.2        # Fraction of Season from Initial to Vegetative
            f2 = 0.5        # Fraction of Season from Initial to Reproductive
            f3 = 0.75       # Fraction of Season from Initial to Ripening
            EoS = 1.0       # Fraction of Season at End

            kc_ini =        # Kc at Initial Stage
            kc_max =        # Kc at Reproductive Stage
            kc_EoS =        # Kc at End of Season

        """
        if not day_of_season >= 0:
            raise ValueError ("day_of_season must be >= 0")
        if day_of_season <= t_seas*f1:
            return kc_ini
        elif day_of_season < t_seas*f2:
            return ((kc_max-kc_ini)/(f2*t_seas-f1*t_seas))*(day_of_season-f1*t_seas)+kc_ini
        elif day_of_season <= t_seas*f3:
            return kc_max
        elif day_of_season < t_seas*EoS:
            return kc_ini+((day_of_season-EoS*t_seas)/(f3*t_seas-EoS*t_seas))*kc_EoS+kc_ini
        else:
            return kc_EoS
            #return self.kc_max # previously this was the last line of the function



    
    def calc_T_max(self, t):
        """ Calculates max Transpiration variable.
        
        Usage: calc_T_max(t)

            t = day of season

        T = kc(t) * T_max
            
        """
        return self.calc_kc(t) * self.T_max

    def calc_LAI(self, day_of_season, p=1):
        """ Returns a Leaf Area Index (LAI) variable. LAI comes
            from function of kc. Currently based on linear relationship 
            between kc and LAI (assumption).
        
        Usage: calc_LAI(t, p=1)

        LAI = (LAI_max/kc_max)^p * kc(t),

        where kc varies through the season according to calc_kc(t) 

        Note: p=1 assumes a linear relationship between LAI and kc

        """
        return pow((self.LAI_max/self.kc_max),p) * self.calc_kc(day_of_season)

    def calc_T(self, s, t):
        """ Calculates Transpiration variable as a stepwise
            linear function.
        
        Usage: calc_T(s, t)
        
            s = relative soil moisture [0-1]
            t = day of season

        """
        if not s <= 1:
            raise ValueError("Time stemp {t}: s ({s}) must be <= 1".format(
                t=t,
                s=s))
        if s>=self.s_star:
            return self.calc_T_max(t)
        elif s>=self.sw:
            return (s-self.sw)/(self.s_star-self.sw)*self.calc_T_max(t)
        else:
            return 0

