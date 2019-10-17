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


class TestModel(unittest.TestCase):

    def setUp(self):
        climate = Climate()
        soil = Soil('sand')
        crop = Crop(soil=soil)
        self.model = CropModel(crop=crop, soil=soil, climate=climate)

    def test_cropModel(self):
        assert isinstance(self.model.R, np.ndarray), "R value is an ndarray"