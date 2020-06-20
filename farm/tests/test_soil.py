"""
Name:           test_soil.py
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

class TestSoil(unittest.TestCase):

    def setUp(self):
        climate = Climate()
        soil_ls = Soil('loamy sand')
        soil_l = Soil('loam')
        crop_ls = Crop(soil=soil_ls)
        crop_l = Crop(soil=soil_l)
        self.model_ls = CropModel(crop=crop_ls,soil=soil_ls,climate=climate)
        self.model_l = CropModel(crop=crop_l,soil=soil_l,climate=climate)

    # 1. Do the hygroscopic values correspond with expected values from Liao et al. 2001b, Figure 6
    #    ca. 0.07 for loamy sand, and 0.19 for loam?
    def test_Liao_ls(self):
        assert 0.02 <= self.model_ls.soil.sh <= 0.12, "Error: sh for loamy sand is too high; it should be around 0.07"

    def test_Liao_l(self):
        assert 0.13 <= self.model_l.soil.sh <= 0.25, "Error: sh for loam is too high; it should be around 0.19"