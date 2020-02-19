#%% SCRIPT TO RUN THE MODEL
#
# Import matplotlib to generate figures
import matplotlib.pyplot as plt

# This is an example script that runs
# the ecohydrology model.
from farm import Climate
from farm import Soil
from farm import Crop
from farm import CropModel

#%% Initialize Climate
climate = Climate()

#%% Initialize Soil
soil = Soil(texture='loam')

#%% Initialize Plant
crop = Crop(soil=soil)

# Finally, now that the plant is defined, we need
# to set the soil.nZr property:
soil.set_nZr(crop)

#%% Initialize the Model
model = CropModel(crop=crop, climate=climate, soil=soil)

#%% Run the Model
model.run()

#%% Assign Model Output
output = model.output()

#%% Generate figures
output['ET'].plot()