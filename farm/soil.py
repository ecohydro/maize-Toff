#%% Soil Class Definition
from math import exp
from math import log as ln # This is kind of janky! Could use numpy instead

"""
This is the "example" module.

The example module supplies one function, factorial().  For example,

>>> Soil('Sand')
<__main__.Soil object at 0x7fe13f96f0b8>

# TODO: Maybe move soil parameters out of this class. 
"""

#%% Set parameters related to soils and siginificant digits
rho = 1000.0    # density of water in kg/m^3
g = 9.8         # acceleration of gravity in m/s^2
PRECISION = 2   # Number of decimal places of precision in calculations (default is 2)

#%% DATA FROM CLAPP AND HORBERGER (C&H) 1978, Table 2:

field_capacity = -33 / 1000 # Field capacity in MPa.

soils = {
    'sand':{
        'b': 4.05,
        'Psi_S_cm': 12.1,   # saturated water tension, cm
        'Psi_l_cm': 4.66,   # leakage water tension, cm
        'n': 0.395,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 1.056,        # saturated hydraulic conductivity, cm/min
        'S': 1.52           # sorptivity, cm/min^1/2    
    },
    'loamy sand':{
        'b': 4.38,
        'Psi_S_cm': 9.0,    # saturated water tension, cm
        'Psi_l_cm': 2.38,   # leakage water tension, cm
        'n': 0.410,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.938,        # saturated hydraulic conductivity, cm/min
        'S': 1.04           # sorptivity, cm/min^1/2  
    },
    'sandy loam':{
        'b': 4.90,
        'Psi_S_cm': 21.8,   # saturated water tension, cm
        'Psi_l_cm': 9.52,   # leakage water tension, cm
        'n': 0.435,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.208,        # saturated hydraulic conductivity, cm/min
        'S': 1.03           # sorptivity, cm/min^1/2  
    },
    'silt loam':{
        'b': 5.30,
        'Psi_S_cm': 78.6,   # saturated water tension, cm
        'Psi_l_cm': 75.3,   # leakage water tension, cm
        'n': 0.485,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0432,       # saturated hydraulic conductivity, cm/min
        'S': 1.26           # sorptivity, cm/min^1/2  
    },
    'loam':{
        'b': 5.39,
        'Psi_S_cm': 47.8,   # saturated water tension, cm
        'Psi_l_cm': 20.0,   # leakage water tension, cm
        'n': 0.451,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0417,       # saturated hydraulic conductivity, cm/min
        'S': 0.693          # sorptivity, cm/min^1/2  
    },
    'sandy clay loam':{
        'b': 7.12,
        'Psi_S_cm': 29.9,   # saturated water tension, cm
        'Psi_l_cm': 11.7,   # leakage water tension, cm
        'n': 0.420,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0378,       # saturated hydraulic conductivity, cm/min
        'S': 0.488          # sorptivity, cm/min^1/2  
    },
    'silty clay loam':{
        'b': 7.75,
        'Psi_S_cm': 35.6,   # saturated water tension, cm
        'Psi_l_cm': 19.7,   # leakage water tension, cm
        'n': 0.477,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0102,       # saturated hydraulic conductivity, cm/min
        'S': 0.310          # sorptivity, cm/min^1/2  
    },
    'clay loam':{
        'b': 8.52,
        'Psi_S_cm': 63.0,   # saturated water tension, cm
        'Psi_l_cm': 48.1,   # leakage water tension, cm
        'n': 0.476,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0147,       # saturated hydraulic conductivity, cm/min
        'S': 0.537          # sorptivity, cm/min^1/2  
    },
    'sandy clay':{
        'b': 10.4,
        'Psi_S_cm': 15.3,   # saturated water tension, cm
        'Psi_l_cm': 8.18,   # leakage water tension, cm
        'n': 0.426,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0130,       # saturated hydraulic conductivity, cm/min
        'S': 0.223          # sorptivity, cm/min^1/2  
    },
    'silty clay':{
        'b': 10.4,
        'Psi_S_cm': 49.0,   # saturated water tension, cm
        'Psi_l_cm': 23.0,   # leakage water tension, cm
        'n': 0.492,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0062,       # saturated hydraulic conductivity, cm/min
        'S': 0.242          # sorptivity, cm/min^1/2  
    },
    'clay':{
        'b': 11.4,
        'Psi_S_cm': 40.5,   # saturated water tension, cm
        'Psi_l_cm': 24.3,   # leakage water tension, cm
        'n': 0.482,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0077,       # saturated hydraulic conductivity, cm/min
        'S': 0.268          # sorptivity, cm/min^1/2  
    }
}

#%% Soil CLass Definition

# pylint: disable=maybe-no-member
class Soil():
    """ Defines a soil object based on either passed parameters
    or a soil texture class corresponding to the textures defined in 
    Clapp & Hornberger (C&H), 1978, Table 2.

    Usage: Soil(texture, params)
        Notes: 
            If texture is not provided, params must be included.
            If texture is provided, params is ignored.
            Capitilization in texture classes is ignored.
        
        texture = texture name from Clapp & Hornberger, 1978, Table 2.
            Valid options are:
            ["Sand", "Loamy Sand", "Sandy Loam", "Silt Loam", "Loam",
             "Sandy Clay Loam", Silty Clay Loam", "Clay Loam", "Sandy Clay",
             "Silty Clay", "Clay"]
        
        params = dictionary containing values for the soil parameters:
            params = {
                'b': 11.4,
                'Psi_S': 40.5,  # saturated water tension, cm
                'Psi_l': 24.3,  # leakage water tension, cm
                'n': 0.482,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
                'Ks': 0.0077,   # saturated hydraulic conductivity, cm/min
                'S': 0.268      # sorptivity, cm/min^1/2  
            }
    
    Note: In C&H 1978, soil water retention relationships were defined according to _tensions_. 
    These tensions are specified as lengths, and are always _positive_ 
    (tension, like depth has an implied relationship to zero).

    To convert a tension, Psi_cm, (positive quantity of length) into a water potential, Psi_Pa, 
    (negative measure of energy density per unit volume, or Pa), you do the following:

    Psi_Pa = -1 * Psi_cm * rho * g

    This conversion is done during initiation of the soil class.

    """
    def __init__(self, texture=None, params=None):
        """ Initializes a soil object.

        The init function requires _either_ a soil texture or a params dictionary
        (see class description)
        
        """
        self._valid_params = set(['b', 'Psi_S_cm', 'Psi_l', 'n', 'Ks', 'S'])
        self._required_params = set(['b', 'Psi_S_cm', 'n', 'Ks'])
        
        # Set required attributes to None:
        [setattr(self, attr, None) for attr in self._required_params]
        
        if texture: # If this class is instanced with a specific USDA soil texture.
            texture = texture.lower() # Force the soil texture category to lower case
            # Assign texture parameters based on the appropriate soil class:
            for attr, val in soils[texture].items():
                setattr(self, attr, val)
        elif params: # If the class is instanced with a set of soil parameters
            for attr, val in params.items():
                # Only include valid soil parameters
                if attr in self._valid_params:  
                    setattr(self, attr, val)
            # Check that all required parameters have been set
            if not self._required_params.issubset(self.__dict__.keys()):
                missing = self._required_params.difference(self.__dict__.keys())
                raise AttributeError("Missing required parameters, {list}".format(list=missing))
        else: 
            raise AttributeError("Must pass either a soil texture or dict of parameters")
        
        # Set Psi_S (MPa) from Psi_S_cm (cm). Assumes that Psi_S_cm is positive (as it should be!)
        self.Psi_S_MPa = -1 * self.Psi_S_cm / 100 * rho * g / 1E6 
        
        # This version of sfc calculation comes from Laio et al. 2001b. Specifically, cf. the discussion
        # on p.714, and equation 15. 
        self.sfc = pow(0.05/(self.Ks*60*24),1/(2*self.b+3))  # Convert Ks in cm/day 
        self.sfc = self.s(psi=field_capacity)
        # Make sure that field capacity is always lower than soil porosity.
        if self.sfc > 1:
            raise ValueError("soil field capacity, {sfc} is larger than 1".format(
                sfc=self.sfc
            ))
        # Hygroscopic point is when soil is so dry no further evaporation will occur.
        self.sh = self.s(self.theta(-12))               # Hygroscopic point in relative soil moisture [0-1] # TODO: We also set this to -100 for testing
        self.nZr = None                                 # TODO: Hygroscopic point is a wonky parameter stuck in the middle code.. consider setting elsewhere

    def _check_nZr(self):
        error = "Error: Calculation depends on value of self.nZr before calling self.set_nZr"
        if not self.nZr:
            raise AttributeError(error)

    def _check_theta(self, theta):
        error = "theta, {theta}, must be be in the interval (0,{n}]".format(
                theta=theta, n=self.n)
        if theta > self.n or theta < 0:
            raise ValueError(error)

    def psi(self, theta):
        """ Return water potential in Pa based 
        on volumetric soil water content in m^3/m^3

        Note: Assumes that Psi is a water potential, and therefore Psi < 0 for unsaturated soils!

        Usage: psi(theta):
        
            theta = soil water content [m^3/m^3]
        
        """
        self._check_theta(theta)
        s = self.s(theta=theta)          
        return round(self.Psi_S_MPa * pow(s,-self.b),PRECISION)
    
    def theta(self,psi):
        """ Return a volumetric water content in m^3/m^3 
        based on a given water potential (MPa)

        Note: 
        Usage: theta(psi):

            psi = soil water potential [MPa]
        
        """
        if psi > 0:
            raise ValueError("psi, {psi}, must be less than or equal to zero.".format(psi=psi))
        # Ensure result is rounded to correct precision and that we do not exceed porosity
        return min([round((self.n * pow(psi/self.Psi_S_MPa, 1/-self.b)),PRECISION), self.n]) 


    def s(self,theta=None,psi=None):
        """ Return a relative soil moisture value, s [0-1]
        given a volumetric water content [m^3/m^3] or a 
        water potential [MPa]

        Usage: s(theta):

            theta = volumetric water content [m^3/m^3]
            psi = water potential [MPa]

        Note: theta must be in the interval 0-n (porosity)
        Note: psi must be negative
        Note: Function must be called with either theta or psi, but not both.
        
        """
        if theta and psi:
            raise ValueError(
            "Both theta ({theta}) and psi {psi} values provided only one argument allowed".format(
                theta=theta,
                psi=psi
            ))
        if psi:
            theta = self.theta(psi)
        self._check_theta(theta)
        try:
            return round(theta/self.n, PRECISION)
        except:
            raise ValueError("Either theta or psi must be provided as an argument.")

    def set_nZr(self,plant):
        """ Sets the nZr for this soil in order to 
        determine fluxes in mm/day rather than relative
        soil moisture

        Usage: set_nZr(plant)

            plant = plant object with plant.Zr value set.
        
        Returns:

            nZr = n * Zr

            Also sets internal soil property nZr according to:

                self.nZr = self.n * plant.Zr
        """
        self.nZr = self.n * plant.Zr 
        return self.nZr 

    
    def calc_Q(self,s,units='mm/day'):
        """ Determines runoff as a function of relative soil moisture

        Usage: 

            calc_Q(s,units)

            s = relative soil moisture [0-1]
            units = units to return leakage in
                options are 'mm/day' (default). 
                Otherwise, returns in [0-1] relative soil 
                moisture

        Returns:

            Q = runoff [mm/day] or [0-1]
        
        """

        # Saturation excess runoff occurs when 
        # relative soil moisture exceeds 1.
        Q = 0
        if s > 1:
            Q = s - 1
        if units == 'mm/day':
            self._check_nZr()           
            return Q * self.nZr
        else:
            return Q

    def calc_L(self,s,units='mm/day'):
        """ Calculates leakage loss as a function of relative soil moisture


        Usage: calc_L(s,units)

            s = relative soil moisture [0-1]
            units = units to return leakage in
                options are 'mm/day' (default). 
                Otherwise, returns in [0-1] relative soil 
                moisture
        Returns:

            L(s) [mm/day] if units='mm/day'
            else returns [0-1]
        
        Notes:
            v1. All soils are assumed to drain to field capacity each day.
            v2. TODO: We determine the daily Leakage using the K(s) function
                and integrating over a day.

        """
        L = 0
        if s > self.sfc:
            L =  min(1, s)-self.sfc
        if units == 'mm/day':
            self._check_nZr()
            return L * self.nZr
        else:
            return L


    def calc_LQ(self, s, t, plant, units='mm/day'):
        """ Calculates the time evolution of soil moisture
        # previously in arguments (not needed?): beta = None, eta = None, mu = None, 
        # Is it necessary to give crop to plant?

        Usage: calc_LQ(s, plant=crop, units)

            # TODO: All of this needs to be fixed
            # TODO: Rename t as day_of_season? day_of_season
            
            # Calculate B for b related to soil type.
            beta = 2 * b + 4 

            # Define eta and mu
            eta_w = E_w / nZr
            eta = E_max / nZr
            mu = K_s / nZr * ( e^beta (1 - s_fc) - )

            # Calculate time variables
            t_sfc = 1 / ( beta * (mu - eta) ) * (
                    beta * (s_fc - s_0) + 
                    ln( eta - mu + mu * e^( beta(s_0 - s_fc) ) / eta ) )

            t_s_star = (s_fc - s_star) / eta + t_sfc

            t_sw = (s_star - s_w) / (eta - eta_w) * ln(eta / eta_w) + t_s_star

            s(t) = 
                0 <= t < t_sfc ... next part

        
        Returns:

            s(t) = soil moisture [mm/day] if units='mm/day' 

        Notes:
            TODO: We have not defined t.
            CHECK: s_0 is the same as s
            CHECK: b is a property of soil above and needs to be calculated for each soil: B

        """
        # #TODO: Start off with S_0 as 0, need to figure this out if that's the case
        s_0 = 0

        # Calculate B for b related to soil type.
        beta = 2 * self.b + 4 


        # TODO: Fix climate variables
        E_w = 3
        E_max = 4

        # TODO: Is s_0 the same as s? If so, change self.s to self.s_0, previously ran into the error:
        # AttributeError: 'Soil' object has no attribute 's_0'
        # Got the same error with self.s0

        # Define eta and mu
        eta_w = E_w / self.nZr # TODO: We do not have a variable for E_w so it is standing in as 3.

        eta = E_max / self.nZr # TODO: E_max comes from climate; how do I get that in here?

        mu =  self.Ks / (self.nZr * ( exp( beta*(1 - self.sfc) ) - 1 ))

        # Calculate time variables
        t_sfc = 1 / ( beta * (mu - eta) ) * (
                    beta * (self.sfc - s_0) + #TODO: I'm not sure if this "s" is correct.
                    ln( eta - mu + mu * exp( beta * (s_0 - self.sfc) ) / eta ) )

        # TODO: s_star is a property of plant...
        t_s_star = (self.sfc - plant.s_star) / eta + t_sfc

        t_sw = (plant.s_star - plant.sw) / (eta - eta_w) * ln(eta / eta_w) + t_s_star

        # Now we can go onto the s(t) calculations...
        if not t >= 0:
            raise ValueError ("t must be >= 0")
       
        if t <= t_sfc:
            return s -  1 / beta * ln( 
                (eta - mu + mu * exp(beta*(s_0 - self.sfc)) * exp( beta * (eta - mu ) * t ) 
                - mu * exp( beta * (s_0 - self.sfc)) 
                )
                / (eta - mu )
                )
        # TODO: Finish the conditions

if __name__ == "__main__":
    import doctest
    doctest.testmod()  