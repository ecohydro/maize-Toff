

rho = 1000  # density of water in kg/m^3
g = 9.8     # acceleration of gravity in m/s^2

#%% DATA FROM CLAPP AND HORBERGER (C&H) 1978, Table 2:

soils = {
    'sand':{
        'b': 4.05,
        'Psi_S': 12.1,  # saturated water tension, cm
        'Psi_l': 4.66,  # leakage water tension, cm
        'n': 0.395,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 1.056,    # saturated hydraulic conductivity, cm/min
        'S': 1.52       # sorptivity, cm/min^1/2    
    },
    'loamy sand':{
        'b': 4.38,
        'Psi_S': 9.0,  # saturated water tension, cm
        'Psi_l': 2.38,  # leakage water tension, cm
        'n': 0.410,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.938,    # saturated hydraulic conductivity, cm/min
        'S': 1.04       # sorptivity, cm/min^1/2  
    },
    'sandy loam':{
        'b': 4.90,
        'Psi_S': 21.8,  # saturated water tension, cm
        'Psi_l': 9.52,  # leakage water tension, cm
        'n': 0.435,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.208,    # saturated hydraulic conductivity, cm/min
        'S': 1.03       # sorptivity, cm/min^1/2  
    },
    'silt loam':{
        'b': 5.30,
        'Psi_S': 78.6,  # saturated water tension, cm
        'Psi_l': 75.3,  # leakage water tension, cm
        'n': 0.485,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0432,    # saturated hydraulic conductivity, cm/min
        'S': 1.26       # sorptivity, cm/min^1/2  
    },
    'loam':{
        'b': 5.39,
        'Psi_S': 47.8,  # saturated water tension, cm
        'Psi_l': 20.0,  # leakage water tension, cm
        'n': 0.451,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0417,    # saturated hydraulic conductivity, cm/min
        'S': 0.693       # sorptivity, cm/min^1/2  
    },
    'sandy clay loam':{
        'b': 7.12,
        'Psi_S': 29.9,  # saturated water tension, cm
        'Psi_l': 11.7,  # leakage water tension, cm
        'n': 0.420,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0378,    # saturated hydraulic conductivity, cm/min
        'S': 0.488       # sorptivity, cm/min^1/2  
    },
    'silty clay loam':{
        'b': 7.75,
        'Psi_S': 35.6,  # saturated water tension, cm
        'Psi_l': 19.7,  # leakage water tension, cm
        'n': 0.477,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0102,    # saturated hydraulic conductivity, cm/min
        'S': 0.310       # sorptivity, cm/min^1/2  
    },
    'clay loam':{
        'b': 8.52,
        'Psi_S': 63.0,  # saturated water tension, cm
        'Psi_l': 48.1,  # leakage water tension, cm
        'n': 0.476,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0147,    # saturated hydraulic conductivity, cm/min
        'S': 0.537       # sorptivity, cm/min^1/2  
    },
    'sandy clay':{
        'b': 10.4,
        'Psi_S': 15.3,  # saturated water tension, cm
        'Psi_l': 8.18,  # leakage water tension, cm
        'n': 0.426,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0130,    # saturated hydraulic conductivity, cm/min
        'S': 0.223       # sorptivity, cm/min^1/2  
    },
    'silty clay':{
        'b': 10.4,
        'Psi_S': 49.0,  # saturated water tension, cm
        'Psi_l': 23.0,  # leakage water tension, cm
        'n': 0.492,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0062,    # saturated hydraulic conductivity, cm/min
        'S': 0.242       # sorptivity, cm/min^1/2  
    },
    'clay':{
        'b': 11.4,
        'Psi_S': 40.5,  # saturated water tension, cm
        'Psi_l': 24.3,  # leakage water tension, cm
        'n': 0.482,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0077,   # saturated hydraulic conductivity, cm/min
        'S': 0.268       # sorptivity, cm/min^1/2  
    }
}

#%% Soil CLass Definition

class Soil():
    """ Defines a soil object based on either passed parameters
    or a soil texture class corresponding to the textures defined in 
    Clapp & Hornberger, 1978, Table 2.

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
            'params':{
                'b': 11.4,
                'Psi_S': 40.5,  # saturated water tension, cm
                'Psi_l': 24.3,  # leakage water tension, cm
                'n': 0.482,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
                'Ks': 0.0077,   # saturated hydraulic conductivity, cm/min
                'S': 0.268       # sorptivity, cm/min^1/2  
            }
    """
    def __init__(self, texture=None, params=None):
        if texture:
            self.b = soils[texture]["b"]
            self.Psi_S = soils[texture]["Psi_S"]
            self.Psi_l = soils[texture]["Psi_l"]
            self.n = soils[texture]["n"]
            self.Ks = soils[texture]["Ks"]
            self.S = soils[texture]["S"]
        elif params:
            self.b = params["b"]
            self.Psi_S = params["Psi_S"]
            self.Psi_l = params["Psi_l"]
            self.n = params["n"]
            self.Ks = params["Ks"]
            self.S = params["S"]

    def psi(self,theta):
        """ Return water potential in Pa based 
        on volumetric soil water content in m^3/m^3

        Usage: psi(theta):
        
            theta = soil water content [m^3/m^3]
        
        """

        return self.n
