"""
Name:           __runmodel__.py
Compatibility:  Python 3.7
Description:    Ecohydro model to track evaporation and transpiration

URL:            https://github.com/ecohydro/maize-Toff

Requires:       pandas; matlotlib; numpy; math

AUTHOR:         Natasha Krell & Kelly Caylor
ORGANIZATION:   U.C. Santa Barbara
Contact:        nkrell@ucsb.edu

"""
    
#%% Import Object Classes
# This is super-weird, in that we are having to reference these from the root dir.
# It will likely break later on, so look out.
from py.climate import Climate
from py.soil import Soil


#%% IMPORT PACKAGES

from math import exp
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



def calc_E(s, E_max, a=2):
    """ Calculates evaporation variable. 
        Based on relationship between evaporation and 
        soil moisture where maximum soil moisture is maximum 
        evapotranspiration. Returns [mm/unit time].
        
    Usage: calc_E(s, E_max, a=2)
        Note: s must be a single-dimension array
        s = from 0 to 1  # Relative soil water storage [mm/mm]
        E_max = some number # Max daily evapotranspiration [mm/day]
        a = 2 # Assume a is >= 0 

    >>> calc_E(1.1,1)
    Traceback (most recent call last):
        ...
    ValueError: s must be <= 1
    >>> calc_E(0.1, 1, -0.5)
    Traceback (most recent call last):
        ...
    ValueError: a must be >= 0
    >>> calc_E(0.5,2)
    0.5

    """
    if not s <= 1:
        raise ValueError("s must be <= 1")
    if not a >= 0:
        raise ValueError("a must be >= 0")
    return E_max * pow(s,a)

def calc_E_max(LAI, k=0.9, ET_max=6):
    """ Calculates maximum (potential) evaporation variable.
    
    Usage: calc_E_max(LAI, k=0.9, ET_max=6)
        LAI = Between 0 and 3 # Leaf Area Index [unitless]
        k = 0.9 # some Constant [unitless]
        ET_max = 6 # Potential ET [mm/day]

    LAI must not be ridiculously large
    >>> calc_E_max(3.1)
    Traceback (most recent call last):
        ...
    ValueError: LAI must be <= 3

    """
    if not LAI <= 3:
        raise ValueError("LAI must be <= 3")

    return ET_max * exp(-k*LAI)

# def calc_LAI(kc, LAI_max=3.5, kc_max=1.2):
#     """ Calculates Leaf Area Index (LAI) variable. LAI comes
#         from function of kc. Currently based on linear relationship 
#         between kc and LAI (assumption).
    
#     Usage: calc_LAI(kc, LAI_max, kc_max)
#         kc = some number # Crop coefficient
#         LAI_max = 3.5 # Max. LAI
#         kc_max = 1.2 # Max. crop coeff.

#     >>> calc_LAI(1)
#     2.916666666666667
    
#     """
#     return (LAI_max/kc_max) * kc

# def calc_T_max(kc, ET_max=6):
#     """ Calculates max. Transpiration variable.
    
#     Usage: calc_T_max(kc, ET_max=6)
#         kc = some number # Crop coefficient
#         ET_max = 6 # max. Evapotranspiration [mm/day]
#     >>> calc_T_max(1)
#     6

#     """
#     return kc * ET_max

# def calc_T(s, T_max, sw=0.3, s_star=0.6):
#     """ Calculates Transpiration variable as a stepwise
#         linear function.
    
#     Usage: calc_T(s, T_max, sw=0.3, s_star=0.6)
#         Note: s must be a single-dimension array
#         s = some value # Relative soil moisture [mm/mm]
#         T_max = defined by a function [mm/unit time]
#         sw = 0.3 # Wilting point [mm]
#         s_star = 0.6 # Plant water satisfication [mm]
    
#     >>> calc_T(1.2,1)
#     Traceback (most recent call last):
#         ...
#     ValueError: s must be <= 1
#     >>> calc_T(0.5,1)
#     0.6666666666666667
#     >>> calc_T(0.2, 1, sw=0.3)
#     0
    
#     """
#     if not s <= 1:
#         raise ValueError("s must be <= 1")
#     if s>=s_star:
#         return T_max
#     elif s>=sw:
#         return (s-sw)/(s_star-sw)*T_max
#     else:
#         return 0

def calc_ET(s,plant=plant,climate=climate):
    """ Calculates Evapotranspiration (finally!)
    
    Usage: calc_ET(s, plant=plant, climate=climate)
        Note: s must be a single-dimension array
        s = user input # Relative soil moisture [mm/mm]
        plant = plant or crop object
        climate = climate object
    
    """
    T = plant.calc_T(s)
    E = climate.calc_E(s)
    return E, T, E+T

# def calc_kc(fraction_of_season=None):
#     """ Add comments """
#     # TODO: Fix this so kc varies throughout season.
#     return 1.0

# def calc_kc(x, mean=60, stand=20, pi=3.14159):
#     """ Meant to calculate kc based on DOY but not
#     yet working.

#     Usage: calc_kc(x, mean=60, stand=20, pi=3.14519)
#         Note: x must be a single-dimension array
#         mean = 60 # average
#         stand = 20 # standard deviation
#         pi = 3.14159 # Pi

#     >>> calc_kc(1)
#     0.000257132046152697
#     """
#     return 1/(stand * np.sqrt(2 * np.pi)) * np.exp( - (x - mean)**2 / (2 * stand**2) )


#%% MAIN MODEL CODE
# 
# dS/dt = R - E - T 
# TODO Add Model equations
#
######################################


#%% Define climate
climate = Climate(
    alpha=10,       # average storm depth [mm]
    lambda_r=0.3,   # storm frequency [day^-1]
    t_seas=120,     # length of season [days]
    ET_max=6.5      # maximum daily ET [mm]
    )



#%% Define the CropModel Class
class CropModel():
    """ Defines an ecohydrological model class for use with crops.


    """
    def __init__(self, *args, **kwargs):
        from numpy import zeros
        self.soil = kwargs.pop('soil')
        self.crop = kwargs.pop('crop')
        self.climate = kwargs.pop('climate')
        self.n_days = len(climate.rainfall)

        # Set maximum soil water content (mm):
        self.nZr = self.soil.n * self.crop.Zr # [mm]  

        # Pre-allocate arrays
        self.R = climate.rainfall
        self.s = zeros(self.n_days)
        self.ET = zeros(self.n_days)
        self.E = zeros(self.n_days)
        self.T = zeros(self.n_days)
        self.L = zeros(self.n_days)
        self.dsdt = zeros(self.n_days)
        
        self.LAI = zeros(self.n_days)
        self.ET_max = zeros(self.n_days)
        self.T_max = zeros(self.n_days)
        self.kc = zeros(self.n_days)

        # Set initial conditions:
        self.s[0] = 0.3     # relative soil moisture, [0-1]
        
    def run(self):
        for t in range(self.n_days):
            try:
                # 0. Update the crop coefficient
                # TODO: Edit Crop class to make this dynamic.
                self.kc[t] = self.crop.calc_kc(t)
                self.LAI[t] = self.crop.calc_LAI(t)

                # 1. Calculate ET terms
                self.T[t] = self.crop.calc_T(self.s[t],t)   # mm/day
                self.E[t] = self.climate.calc_E(
                    self.s[t],
                    t,
                    plant=self.crop,
                    soil=self.soil) # mm/day
                self.ET[t] = self.T[t] + self.E[t]
                
                # 2. Update Soil Moisture Water Balance
                self.dsdt[t] = self.R[t] - self.ET[t]           # mm/day
                self.s[t+1] = self.s[t] + self.dsdt[t]/self.nZr # [0-1]
                
                # 3. Force Soil Moisture Water Balance to Limits
                if self.s[t+1] > 1:
                    self.dsdt[t] = self.dsdt[t] - (self.s[t+1] - 1)
                    self.L[t]  = (self.s[t+1] - 1)*self.nZr # [mm/day]
                    self.s[t+1] = 1
                if self.s[t+1] < 0:
                    self.s[t+1] = 0
                
                # 4. Determine leakage loss:
                if self.s[t+1] > self.soil.sfc:
                    # Update dsdt for leakage loss
                    self.dsdt[t]  = self.dsdt[t] - (self.s[t+1] - self.soil.sfc)
                    # Calculate leakage in units of mm/day
                    self.L[t] = (self.s[t+1] - self.soil.sfc)*self.nZr # [mm/day]
                    self.s[t+1] = self.soil.sfc
            except:
                break
    
    def output(self):
        from pandas import DataFrame
        return DataFrame({ 
            'kc':self.kc,
            'LAI':self.LAI,
            'R':self.R,
            's':self.s,
            'E':self.E,
            'T':self.T,
            'L':self.L,
            'dsdt':self.dsdt,

        })



                   

#%% OLD MODEL CODE

def runmodel(n=0.4, Zr=500, ET_max=6.5, R):
    # INITIALIZE PARAMETERS FOR MODEL
    n = 0.4  # Porosity, [m3/m3]
    Zr = 500 # Rooting depth [mm]
    S_max = n*Zr    # Max soil water storage [mm]
    S_fc = 0.75 * S_max # Field capacity of soil water [mm]
    ET_max = 6.5 # Daily reference evapotranspiration [mm]
    
    # STEP 2. INITIALIZE WATER BALANCE TERMS
    S = np.zeros(len(R)) # Pre-allocate the array for soil moisture
    ET = np.zeros(len(R)) # Pre-allocate the array for evapotranspiration
    L = np.zeros(len(R)) # Pre-allocate the array for leakage loss from soil
    S[0] = 30  # Initial soil water storage [mm]

    dSdt = np.zeros(len(R))# Pre-allocate the array for dSdt

    for t in range(len(R)-1): # Range starts with zero so ignore the last day.
        
        # Update the crop coefficient
        crop.calc_kc(t)

        # Calculate ET
        T[t] = crop.calc_T(s)
        E[t] = climate.calc_E(s,) 
        
        # Update Soil Moisture Water Balance
        dSdt[t] = R[t] - ET[t] # We will handle leakage  (L[t]) later.
        S[t+1] = S[t] + dSdt[t]

        # Force Soil Moisture Water Balance to Limits
        # Soil Moisture must follow: 0 < S <= S_max
        if S[t+1] > S_max:
            L[t] = S[t+1] - S_max # Leakage is excess soil moisture.
            # TODO: Include leakage losses between S_max and S_fc (field capacity)
            S[t+1] = S_max
        if S[t+1] < 0:
            S[t+1] = 0

        # Determine Leakage losses:
        # Assume that the soil always drains to field capacity within a single day.
        # This means that L[t] = S[t+1] - S_fc 
        if S[t+1] > S_fc:
            L[t] = L[t] + S[t+1] - S_fc

    # Package up results and return them all as a dataframe
    result =  pd.DataFrame({
        'Soil Moisture': S,
        'ET': ET,
        'Lakeage Loss': L,
       'Water Balance': dSdt})

    rf_result = pd.DataFrame(R.T, columns=['Value'])
    #return result
    return rf_result.head()

    # TODO: Specify whether to return result or rf_result

    #return R, S, ET, L, dSdt

# if __name__ == "__main__":
    #import doctest
    #doctest.testmod()

#%% MAKE FIGURES
