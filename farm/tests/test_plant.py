"""
Name:           test_plant.py
Compatibility:  Python 3.7
Description:    Integration test of the model

URL:            https://github.com/ecohydro/maize-Toff
"""
import unittest
import numpy as np

from farm import Climate
from farm import Soil
from farm import Crop
from farm import CropModel
from farm.functions import *

class TestPlant(unittest.TestCase):
    def setUp(self):
        climate = Climate()
        soil_l = Soil('loam')
        crop_l = Crop(soil=soil_l)
        self.model = CropModel(crop=crop_l,soil=soil_l,climate=climate)
        

	# Write tests
	# 1. Is the average static stress calculation correct?
    def test_mstr_memb(self):
        self.model.run()
        o = self.model.output()

        mstr_memb, dstr_memb, yield_kg_ha = self.model.crop.calc_dstress(o.s, o.stress)

        # Calculate average stress
        QPAR_CROP = 2
        
        sstr_memb = ((self.model.crop.s_star - o.s) / (self.model.crop.s_star - self.model.crop.sw))  # dim
        sstr_memb = sstr_memb[sstr_memb > 0.]
        sstr_memb[sstr_memb > 1.] = 1.

		# Subset model output to just the growing season
        start = 21
        end = start + self.model.crop.lgp
        sstr_memb_subset = sstr_memb[start:end]

		# Get the average static stress
        if len(sstr_memb) > 0:
            _mstr_memb = np.mean(sstr_memb_subset**QPAR_CROP)  # dim
        else:
            _mstr_memb = 0.  # dim

        assert mstr_memb == _mstr_memb, "Error: Average static stress calculation in model should be the same as _mstr_memb calculated in this unittest."

if __name__ == '__main__':
   unittest.main()

