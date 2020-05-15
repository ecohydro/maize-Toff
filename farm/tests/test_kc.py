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

    # write tests
    def inspect_kc(self):
        assert max(model.kc) == 1.2, "Error: The max value for kc should be 1.2."
