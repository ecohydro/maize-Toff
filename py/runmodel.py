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


#%% IMPORT PACKAGES

from math import exp
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#%% DEFINE FUNCTIONS

def calc_R(alpha, lambda_r, t_seas, n_seasons=1):
    """
    Calculates rainfall variable. 
    Assume an exponential distribution of storm depth. Returns
    amounts and number of rainy days [mm, days].
    
    Usage: calc_R(alpha, lambda_r, t_seas, n_seasons=1)
        alpha = any number # Definition
        lambda_r = any number # Definition
        t_seas = some number greater than 0
        n_seasons = some number greater than 0 # number of seasons

    >>> calc_R()
    Traceback (most recent call last):
        ...
    TypeError: calc_R() missing 3 required positional arguments: 'alpha', 'lambda_r', and 't_seas'
    >>> calc_R(10,0,120,500)
    Traceback (most recent call last):
      ...
    TypeError: object of type 'int' has no len()
    
    """
    amounts = np.random.exponential(scale=alpha, size=[n_seasons, t_seas])
    
    # print('{lam} is not 1 and length {lam_len} is not {t_seas}'.format(
    #    lam=lambda_r,
    #    lam_len=len(lambda_r),
    #    t_seas=t_seas))

    if not isinstance(lambda_r, float):
        if len(lambda_r) != t_seas:
            raise ValueError("lamda_r values should be a constant, or have length of t_seas")
    
    rain_days = (np.random.uniform(low=0, high=1, size=[n_seasons, t_seas]) <= lambda_r).astype(int)
    return amounts * rain_days


#def valid_s(s):


#   return True
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

def calc_LAI(kc, LAI_max=3.5, kc_max=1.2):
    """ Calculates Leaf Area Index (LAI) variable. LAI comes
        from function of kc. Currently based on linear relationship 
        between kc and LAI (assumption).
    
    Usage: calc_LAI(kc, LAI_max, kc_max)
        kc = some number # Crop coefficient
        LAI_max = 3.5 # Max. LAI
        kc_max = 1.2 # Max. crop coeff.

    >>> calc_LAI(1)
    2.916666666666667
    
    """
    return (LAI_max/kc_max) * kc

def calc_T_max(kc, ET_max=6):
    """ Calculates max. Transpiration variable.
    
    Usage: calc_T_max(kc, ET_max=6)
        kc = some number # Crop coefficient
        ET_max = 6 # max. Evapotranspiration [mm/day]
    >>> calc_T_max(1)
    6

    """
    return kc * ET_max

def calc_T(s, T_max, sw=0.3, s_star=0.6):
    """ Calculates Transpiration variable as a stepwise
        linear function.
    
    Usage: calc_T(s, T_max, sw=0.3, s_star=0.6)
        Note: s must be a single-dimension array
        s = some value # Relative soil moisture [mm/mm]
        T_max = defined by a function [mm/unit time]
        sw = 0.3 # Wilting point [mm]
        s_star = 0.6 # Plant water satisfication [mm]
    
    >>> calc_T(1.2,1)
    Traceback (most recent call last):
        ...
    ValueError: s must be <= 1
    >>> calc_T(0.5,1)
    0.6666666666666667
    >>> calc_T(0.2, 1, sw=0.3)
    0
    
    """
    if not s <= 1:
        raise ValueError("s must be <= 1")
    if s>=s_star:
        return T_max
    elif s>=sw:
        return (s-sw)/(s_star-sw)*T_max
    else:
        return 0

def calc_ET(s, kc, ET_max=6, sw=0.3, s_star=0.6):
    """ Calculates Evapotranspiration (finally!)
    
    Usage: calc_ET(s, kc, ET_max=6, sw=0.3, s_star=0.6)
        Note: s must be a single-dimension array
        s = user input # Relative soil moisture [mm/mm]
        kc = user input # Crop coefficient [Unitless]
        ET_max = 6 # Max. Evapotranspiration [mm/day]
        sw = 0.3 # Soil wilting point [mm]
        s_star = 0.6 # Plant water satisfaction [mm]

    >>> calc_ET(1.2,1)
    Traceback (most recent call last):
        ...
    ValueError: s must be <= 1
    >>> calc_ET(0.5,1)
    (0.10865963555137714, 4.0, 4.108659635551377)
    >>> calc_ET(0.2, 1, sw=0.3)
    (0.017385541688220346, 0, 0.017385541688220346)
    
    """
    LAI = calc_LAI(kc)
    T_max = calc_T_max(kc, ET_max)
    E_max = calc_E_max(LAI)
    T = calc_T(s, T_max, sw, s_star)
    E = calc_E(s, E_max)
    return E, T, E+T

def calc_kc(fraction_of_season=None):
    """ Add comments """
    # TODO: Fix this so kc varies throughout season.
    return 1.0

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

    


def runmodel(n=0.4, Zr=500, ET_max=6.5):
    # INITIALIZE PARAMETERS FOR MODEL
    n = 0.4  # Porosity, [m3/m3]
    Zr = 500 # Rooting depth [mm]
    S_max = n*Zr    # Max soil water storage [mm]
    S_fc = 0.75 * S_max # Field capacity of soil water [mm]
    ET_max = 6.5 # Daily reference evapotranspiration [mm]

    alpha = 10      # average storm depth [mm]
    lambda_r = 0.3  # daily probability of rainfall [day^-1]
    t_seas = 120    # length of season [days]

    # STEP 1. SIMULATE DAILY RAINFALL
    R = calc_R(alpha, lambda_r, t_seas) # Implicit that this is one season

    # STEP 2. INITIALIZE WATER BALANCE TERMS
    S = np.zeros(len(R)) # Pre-allocate the array for soil moisture
    ET = np.zeros(len(R)) # Pre-allocate the array for evapotranspiration
    L = np.zeros(len(R)) # Pre-allocate the array for leakage loss from soil
    S[0] = 30  # Initial soil water storage [mm]

    dSdt = np.zeros(len(R))# Pre-allocate the array for dSdt

    for t in range(len(R)-1): # Range starts with zero so ignore the last day.
        # Calculate ET
        kc = calc_kc(fraction_of_season=t/t_seas) # Returns a constant for now. 
        ET[t] = calc_ET(S[t], kc, ET_max=ET_max) # Daily ET [mm/d]

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
