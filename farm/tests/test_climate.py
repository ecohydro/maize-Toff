"""
Name:           test_climate.py
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


class TestClimate(unittest.TestCase):
    