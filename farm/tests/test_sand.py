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


class TestSand(unittest.TestCase):

    def setUp(self):
        climate = Climate()
        soil = Soil('sand')
        crop = Crop(soil=soil)
        self.model = CropModel(crop=crop, soil=soil, climate=climate)

    # write tests
    def test_sh(self):
        assert self.model.soil.sh <= self.model.crop.sw, "Should be true"

    def test_star(self):
        assert self.model.crop.s_star <= self.model.soil.sfc, "Error: s_star is too high: it should be less than field capacity"