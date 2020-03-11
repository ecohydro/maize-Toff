#%%
class Plant():
    """ Defines a plant class


    """
    def __init__(self,
        Zr=500,             # Rooting depth [mm]
        T_max=4,            # Maximum transpiration [mm/day]
        sw_MPa = -1.5,      # wilting point of plant in water potential [MPa]
        s_star_MPa = -0.05, # water potential of maximum transpiration [MPa], 
        soil=None          # a soil in which this plant will grow
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
    def __init__(self, kc_max=1.2, LAI_max=3.0, T_max=4, lgp = 180, f1 = 0.2, f2 = 0.5, f3 = 0.75, 
                 EoS = 1.0, kc_ini = 0.30, kc_EoS = 0.6,*args,**kwargs):

        self.kc_max = kc_max      # Maximum crop coefficient; Kc at Reproductive Stage [0-1]
        self.LAI_max = LAI_max    # Maximum crop leaf area index [m^2/m^2]
        self.T_max = T_max        # Maximum crop water use [mm/day]
        self.lgp = lgp            # Length of growing period [days]
        self.f1 = f1              # Fraction of Season from Initial to Vegetative
        self.f2 = f2              # Fraction of Season from Initial to Reproductive
        self.f3 = f3              # Fraction of Season from Initial to Ripening
        self.EoS = EoS            # Fraction of Season at End
        self.kc_ini = kc_ini      # Kc at Initial Stage
        self.kc_EoS = kc_EoS      # Kc at End of Season
        super(Crop, self).__init__(*args, **kwargs)

    def calc_kc(self, day_of_season=0):
        """ Calculates crop coefficient that varies throughout the season 
        
        Usage: calc_kc(self, day_of_season)
            
            day_of_season = 0 # Start date [day]

        Note: t must be a single-dimension array. day_of_season is when the plant starts
        to grow and is hard coded to zero because python is zero-indexed.
        Kc values come from Table 11. Lengths of crop development stages for 
        maize (grain) in East Africa (altitude), Chapter 6 in FAO (1998).

        """
        if not day_of_season >= 0:
            raise ValueError ("day_of_season must be >= 0")
        if day_of_season <= self.lgp*self.f1:
            return self.kc_ini
        elif day_of_season < self.lgp*self.f2:
            return ((self.kc_max-self.kc_ini)/(self.f2*self.lgp-self.f1*self.lgp))*(day_of_season-self.f1*self.lgp)+self.kc_ini
        elif day_of_season <= self.lgp*self.f3:
            return self.kc_max
        elif day_of_season < self.lgp*self.EoS:
            return self.kc_ini+((day_of_season-self.EoS*self.lgp)/(self.f3*self.lgp-self.EoS*self.lgp))*self.kc_EoS+self.kc_ini
        else:
            return self.kc_EoS

    def calc_T_max(self, kc):
        """ Calculates max Transpiration variable.
        
        Usage: calc_T_max(kc)

        T = kc * T_max
            
        """
        return kc * self.T_max

    def _kc_from_LAI(self, LAI, p=1):
        """ Returns a kc variable. kc comes
            from function of LAI. Currently based on linear relationship 
            between kc and LAI (assumption).
        
        Usage: _kc_from_LAI(LAI, p=1)

            kc = (kc_max/LAI_max)^p * LAI.
            
        Note: p=1 assumes a linear relationship between LAI and kc

        """
        return pow((self.kc_max/self.LAI_max),p) * LAI

    def calc_LAI(self, kc, p=1):
        """ Returns a Leaf Area Index (LAI) variable. LAI comes
            from function of kc. Currently based on linear relationship 
            between kc and LAI (assumption).
        
        Usage: calc_LAI(kc, p=1)

            LAI = (LAI_max/kc_max)^p * kc.

        Note: p=1 assumes a linear relationship between LAI and kc

        """
        return pow((self.LAI_max/self.kc_max),p) * kc
        
    def calc_T(self, s, LAI=None, kc=None):
        """ Calculates Transpiration variable as a stepwise
            linear function. Will use LAI value if both LAI 
            and kc are provided.
        
        Usage: calc_T(s, LAI, kc)
        
            s = relative soil moisture [0-1]
            LAI = leaf area index [m2/m2]
            kc = crop coefficient [-].

        Note: Either LAI or kc must be provided.

        """
        #if not LAI and not kc:
            #raise(ValueError, "Function requires either LAI or kc to be set.")
        if LAI:
            kc = self._kc_from_LAI(LAI)
        if kc:
            if s>=self.s_star:
                return self.calc_T_max(kc)
            elif s>=self.sw:
                return (s-self.sw)/(self.s_star-self.sw)*self.calc_T_max(kc)
            else:
                return 0

    def calc_stress(self, s, q=2):
        """ Calculates static water stress.

        Usage: calc_stress(s, q=2)

            s = relative soil moisture [0-1]
            q = 2

        Note: The value of q changes based on plant species or soil type. 
        See equation 4.13, p.101 in Rodriguez-Iturbe & Porporato (2004)
        """
        if s < self.sw:
            stress = 1
        elif s >= self.s_star:
            stress = 0
        else:
            stress = ((self.s_star - s)/(self.s_star - self.sw))**q
        return stress


