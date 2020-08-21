"""
Name:           test_kc.py
Compatibility:  Python 3.7
Description:    Integration test of the model

URL:            https://github.com/ecohydro/maize-Toff
"""

import unittest
import numpy as np
import matplotlib.pyplot as plt

from farm import Climate
from farm import Soil
from farm import Crop
from farm import CropModel
from farm.functions import *


class TestKc(unittest.TestCase):

    def setUp(self):
        climate = Climate()
        soil = Soil('sand')
        crop = Crop(soil=soil)
        self.model = CropModel(crop=crop, soil=soil, climate=climate)

        #%% Plot Kc
        d = np.arange(121)
        plt.plot(model.kc, '-')
        plt.title('Relationship between Kc and Time of Season')
        plt.xlabel('Time of Season, $\mathit{t}$')
        plt.show()

    def run_model(self, n_sim = 10, burn_in=60, station = 'JACOBSON FARM', texture = 'clay loam', lgp=180, pd_sim=60):
        # Part 1. Set conditions for IC runs
        doy=abs(pd_sim-burn_in)

        # Part 2. Initialize model with a climate, soil and crop
        climate = Climate(station=station)
        soil = Soil(texture=texture)
        crop = Crop(soil=soil, lgp=lgp)
        soil.set_nZr(crop)

        self.model = CropModel(crop=crop, climate=climate, soil=soil)
        model.run()
        o = model.output()

        # Part 3. Get the mean, SD soil moisture and run the simulations to remove IC
        s0_mean, s0_std = average_soil_moisture(model, n_sims=n_sim, doy=doy)
        models = [CropModel(crop=crop, climate=Climate(), soil=soil) for i in np.arange(n_sim)]
        
        # Part 4. Run the actual simulations
        planting_date=pd_sim
        output = [model.run(s0=s0_mean, do_output=True, planting_date=planting_date) for model in models]

        # Part 5. Subset the growing period because model starts 21 days before planting
        start = burn_in
        end = start + lgp
        kc = output[1][start:end]['kc'].tolist() # get kc values for one simulation

        return kc

    # write tests
    # test that crop coefficient does reach the max, 1.2
    def inspect_kc(self):
        assert max(model.kc) == 1.2, "Error: The max value for kc should be 1.2."

    # check that length of kc is the same as the variety 
    def inspect_length(self):
        np.arange(75,180,10)
        result = []

        # This is good enough if it is changing
        for i in varieties:
            # note: bumping up to 200 sims per cultivar type takes several minutes to run
            # whereas 100 sims where varieties = np.arange(70,200,5) takes less than a minute.
            kc = run_model(n_sim = 10, station = 'OL JOGI FARM', texture = 'clay loam', lgp=i, pd_sim=60)
            
            # here we want a list of "trues" indicating that for each variety the length of the kc time series is the 
            # same as the days to maturity
            result.append(i == len(kc))    
            #print(kc)

        # We want to know if all of the kc lengths equal DTM, so the following should be True.    
        res = all(x == True for x in result)
        assert res == True, "Error: the length of the crop coefficient time series, kc is not the same as the days to maturity for the variety."

    # test that the max crop coefficient does occur at a certain point
