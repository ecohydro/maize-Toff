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
        
    def test_climate_init_with_int(self):
       pass
    