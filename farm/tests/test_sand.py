"""
Name:           test_sand.py
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

# TODO: Consider subsuming test_sand into test_soil, and run these tests for various soil types
class TestSand(unittest.TestCase):

    def setUp(self):
        climate = Climate()
        soil = Soil('sand')
        crop = Crop(soil=soil)
        self.model = CropModel(crop=crop,soil=soil,climate=climate)

	# write tests
	# 1. Is the hygroscopic soil moisture less than soil moisture at wilting point?
    def test_sh(self):
        assert self.model.soil.sh <= self.model.crop.sw, "Should be true"


    # 2. Is field capacity less than 1?
    def test_sfc(self):
        assert self.model.soil.sfc <= 1, "Error: field capacity should be less than 1"

    # 3. Is s* less than field capacity?
    def test_star(self):
        assert self.model.crop.s_star <= self.model.soil.sfc, "Error: s_star is too high; it should be less than field capacity"