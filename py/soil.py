

#%% Set parameters related to soils and siginificant digits
rho = 1000      # density of water in kg/m^3
g = 9.8         # acceleration of gravity in m/s^2
PRECISION = 2    # Number of decimal places of precision in calculations (default is 2)

#%% DATA FROM CLAPP AND HORBERGER (C&H) 1978, Table 2:

soils = {
    'sand':{
        'b': 4.05,
        'Psi_S_cm': 12.1,   # saturated water tension, cm
        'Psi_l_cm': 4.66,      # leakage water tension, cm
        'n': 0.395,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 1.056,        # saturated hydraulic conductivity, cm/min
        'S': 1.52           # sorptivity, cm/min^1/2    
    },
    'loamy sand':{
        'b': 4.38,
        'Psi_S_cm': 9.0,    # saturated water tension, cm
        'Psi_l_cm': 2.38,      # leakage water tension, cm
        'n': 0.410,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.938,        # saturated hydraulic conductivity, cm/min
        'S': 1.04           # sorptivity, cm/min^1/2  
    },
    'sandy loam':{
        'b': 4.90,
        'Psi_S_cm': 21.8,   # saturated water tension, cm
        'Psi_l_cm': 9.52,      # leakage water tension, cm
        'n': 0.435,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.208,        # saturated hydraulic conductivity, cm/min
        'S': 1.03           # sorptivity, cm/min^1/2  
    },
    'silt loam':{
        'b': 5.30,
        'Psi_S_cm': 78.6,   # saturated water tension, cm
        'Psi_l_cm': 75.3,      # leakage water tension, cm
        'n': 0.485,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0432,       # saturated hydraulic conductivity, cm/min
        'S': 1.26           # sorptivity, cm/min^1/2  
    },
    'loam':{
        'b': 5.39,
        'Psi_S_cm': 47.8,   # saturated water tension, cm
        'Psi_l_cm': 20.0,      # leakage water tension, cm
        'n': 0.451,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0417,       # saturated hydraulic conductivity, cm/min
        'S': 0.693          # sorptivity, cm/min^1/2  
    },
    'sandy clay loam':{
        'b': 7.12,
        'Psi_S_cm': 29.9,   # saturated water tension, cm
        'Psi_l_cm': 11.7,      # leakage water tension, cm
        'n': 0.420,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0378,       # saturated hydraulic conductivity, cm/min
        'S': 0.488          # sorptivity, cm/min^1/2  
    },
    'silty clay loam':{
        'b': 7.75,
        'Psi_S_cm': 35.6,   # saturated water tension, cm
        'Psi_l_cm': 19.7,      # leakage water tension, cm
        'n': 0.477,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0102,       # saturated hydraulic conductivity, cm/min
        'S': 0.310          # sorptivity, cm/min^1/2  
    },
    'clay loam':{
        'b': 8.52,
        'Psi_S_cm': 63.0,   # saturated water tension, cm
        'Psi_l_cm': 48.1,      # leakage water tension, cm
        'n': 0.476,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0147,       # saturated hydraulic conductivity, cm/min
        'S': 0.537          # sorptivity, cm/min^1/2  
    },
    'sandy clay':{
        'b': 10.4,
        'Psi_S_cm': 15.3,   # saturated water tension, cm
        'Psi_l_cm': 8.18,      # leakage water tension, cm
        'n': 0.426,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0130,       # saturated hydraulic conductivity, cm/min
        'S': 0.223          # sorptivity, cm/min^1/2  
    },
    'silty clay':{
        'b': 10.4,
        'Psi_S_cm': 49.0,   # saturated water tension, cm
        'Psi_l_cm': 23.0,      # leakage water tension, cm
        'n': 0.492,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0062,       # saturated hydraulic conductivity, cm/min
        'S': 0.242          # sorptivity, cm/min^1/2  
    },
    'clay':{
        'b': 11.4,
        'Psi_S_cm': 40.5,   # saturated water tension, cm
        'Psi_l_cm': 24.3,      # leakage water tension, cm
        'n': 0.482,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0077,       # saturated hydraulic conductivity, cm/min
        'S': 0.268          # sorptivity, cm/min^1/2  
    }
}

#%% Soil CLass Definition

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
                'S': 0.268       # sorptivity, cm/min^1/2  
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
        self.Psi_S_MPa = -1 * self.Psi_S_cm * rho * g / 1E6
        self.Psi_L_MPa = -1 * self.Psi_l_cm * rho * g / 1E6
        self.sfc = self.s(self.theta(self.Psi_L_MPa))   # Field capacity in relative soil moisture [0-1]

        # Hygroscopic point is when soil is so dry no further evaporation will occur.
        self.sh = self.s(self.theta(-12))               # Hygroscopic point in relative soil moisture [0-1]
        self.nZr = None

    def _check_nZr(self):
        error = "Error: Calculation depends on value of self.nZr before calling self.set_nZr"
        if not self.nZr:
            raise AttributeError(error)

    def _check_theta(self, theta):
        error = "theta, {theta}, must be be in the interval (0,{n}]".format(
                theta=theta, n=self.n)
        if theta > self.n or theta <= 0:
            raise ValueError(error)

    def psi(self, theta):
        """ Return water potential in Pa based 
        on volumetric soil water content in m^3/m^3

        Note: Assumes that Psi is a water potential, and therefore Psi < 0 for unsaturated soils!

        Usage: psi(theta):
        
            theta = soil water content [m^3/m^3]
        
        """
        self._check_theta(theta)          
        return round(self.Psi_S_MPa * pow(self.n/theta,self.b),PRECISION)
    
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

    def s(self,theta):
        """ Return a relative soil moisture value, s [0-1]
        given a volumetric water content [m^3/m^3]

        Usage: s(theta):

            theta = volumetric water content [m^3/m^3]

        Note: theta must be in the interval 0-n (porosity)
        
        """
        self._check_theta(theta)
        return round(theta/self.n, PRECISION)

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
            Q = 1 - s
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

        """
        L = 0
        if s > self.sfc:
            L =  min(1, s)-self.sfc
        if units == 'mm/day':
            self._check_nZr()
            return L * self.nZr
        else:
            return L



        

