"""
Name:           test_model.py
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


class TestModel(unittest.TestCase):

    def setUp(self):
        climate = Climate()
        soil_l = Soil('loam')
        crop_l = Crop(soil=soil_l)
        self.model = CropModel(crop=crop_l,soil=soil_l,climate=climate)
    
	# Write tests
	# 1. Check that the output is 180 after accounting for burn in
    def test_length_of_output(self):
        self.model.run()
        o = self.model.output()

        start = 60 # Not sure how to get self.burn_in
        end = start + 180 # lgp shouldn't be hard coded in

        lgp = len(o[start:end])

        return lgp

        assert lgp == 180, "Error: Check the length of your dataframe and that burn_in time is accounted for"

if __name__ == '__main__':
   unittest.main()