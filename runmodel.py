#%% SCRIPT TO RUN THE MODEL
#
# This is an example script that runs 
# the ecohydrology model.
from py.climate import Climate
from py.soil import Soil
from py.plant import Crop
from py.model import CropModel

#%% Setup a Parameter Dictionary
params = {
    'climate': {
        'alpha_r':10,       # Average storm depth [mm]
        'lambda_r': 0.25,   # Average storm frequency [day^-1]
        't_seas': 180,      # Length of rainy season [days]
        'ET_max': 6.5       # Maximum ET [mm/day]
    },
    'soil': {
        'texture':'loam' # USDA soil texture
    },
    'crop': {
        'Zr': 500,          # Planting depth [mm]
        'sw_MPa':-1.5,      # Plant wilting point [MPa]
        's_star_MPa':-0.2,  # Water potential for max T
        'kc_max':1.2,       # Maximum crop coefficient
        'LAI_max':3.0,      # Max Leaf Area Index [m2/m2]
        'T_max':4.0         # Max Crop Water Use [mm/day]
    }
}

#%% Initialize Climate
climate = Climate(
    alpha_r=params['climate']['alpha_r'],
    lambda_r=params['climate']['lambda_r'],
    t_seas=params['climate']['t_seas'],
    ET_max=params['climate']['ET_max']
)

#%% Initialize Soil
soil = Soil(texture=params['soil']['texture'])

#%% Initialize Plant
crop = Crop(
    Zr=params['crop']['Zr'],
    sw_MPa=params['crop']['sw_MPa'],
    s_star_MPa=params['crop']['s_star_MPa'],
    kc_max=params['crop']['kc_max'],
    LAI_max=params['crop']['LAI_max'],
    T_max=params['crop']['T_max'],
    soil=soil
)

#%% Initialize the Model
model = CropModel(crop=crop,climate=climate,soil=soil)

#%% Run the Model
model.run()

#%% Assign Model Output
output = model.output()

