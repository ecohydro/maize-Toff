#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2019'

__license__ = 'MIT'
__date__ = 'Thu Feb 27 12:59:28 2020'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''


"""

Name:           filename.py
Compatibility:  Python 3.7.0
Description:    Description of what program does

URL:            https://

Requires:       list of libraries required

Dev ToDo:       None

AUTHOR:         Bryn Morgan
ORGANIZATION:   University of California, Santa Barbara
Contact:        bryn.morgan@geog.ucsb.edu
Copyright:      (c) Bryn Morgan 2019


"""


#%% IMPORTS

# import packages and set working directory
import numpy as np
import matplotlib.pyplot as plt
import os
from math import exp
import pandas as pd

# Change working directory
os.chdir('/Users/brynmorgan/Research/EcohydroModel/maize-Toff/')


# import objects
from farm import Climate
from farm import Soil
from farm import Crop
from farm import CropModel


#%% FUNCTIONS


def annual_sm(simulations,sm0=0.3):
    # generate little progress bar
    import time
    # import progressbar
    # for i in progressbar.progressbar(range(simulations)):
    #     time.sleep(0.02)
    df = pd.DataFrame()
    #result = []
    # run this many times
    for i in range(simulations):
        # Step 1.  Start DOY before the planting season.
        # Pick a day when there isn’t much rain: July 1, Julian Day 182 in non-leap year.
        # Step 2. Run t_sim for 365 so that calendar ends on June 30.
        new_climate = {
            'alpha_r': [9.247170,9.989583,10.586111,11.529771,9.157616,10.351128,10.967919,9.553691,9.254545,8.326368,8.137118,6.693333],
            #[10.0] * 12,
            'lambda_r': [0.051808,0.051502,0.105572,0.272917,0.304435,0.138542,0.174395,0.150202,0.114583,0.209157,0.246237,0.093652],
            'doy_start': 182,
            't_sim': 365,
            'ET_max': 6.5
        }
        climate = Climate(climate_parameters=new_climate)
        soil = Soil('loam')
        crop = Crop(soil=soil)
        soil.set_nZr(crop)
        # Step 3. Start planting date on some reasonable planting date
        # Picked early April Julian Day 100 in regular year.
        model = CropModel(crop=crop,soil=soil,climate=climate, planting_date=100)
        model.run(s0=sm0)
        fin = model.output()
        
        # Return just the soil moisture in a dataframe and append the simulation
        # number to x
        df['x_' + str(i)] = fin.s
        
        # Step 4. Take the soil moisture of June 30 and
        # Re-run the simulation where the soil moisture on July 1st is the value that ended on June 30.
        #print(fin.s.tail(1)) # The argument in tail is the last n lines
        sm0 = fin.s.tail(1)
        # Step 5. Re-run the simulation where the soil moisture on July 1st is the value that ended on June 30.
        # To do this we want to set s0 in model.run as last_s

        # Step 6. Run this lots and lots of times and return a matrix where the
        # row is each simulation and the column is doy

    # Step 7. Remove the first year of data because that was an artificial soil moisture value.
    # So drop the first column (this needs to be out of the loop)
    #df = df.drop(['x_0'], axis=1)
    return df # could transpose so that simulation is row and DOY is column






#%% MAIN
def main():
    
    
#%%


    result = annual_sm(50)
    result
    # Step 8. End goal: For a given day of the year, return the mean and
    # variance of the soil moisture for that/each day. Plot the results.
    df = result
    df['mean'] = df.mean(axis=1)
    df['var'] = df.var(axis=1)
    df

    
    # This is for one climatology and one maize variety (180 days)
    fig = plt.figure()
    plt.plot(df['mean'], 'k-')
    x = np.arange(0, 365, 1) #np.arange(90, 135, 1)
    
    plt.fill_between(x, df['var']+df['mean'], df['mean']-df['var'],facecolor='lightblue') #, facecolor=‘lightblue’
    plt.title('Annual soil moisture for 180 day variety \n Planted on Julian Day 100')
    
    plt.ylabel('Relative soil moisture content, $\mathit{s}$')
    plt.xlabel('Day of season, $\mathit{d}$')
    
    plt.legend(['Mean', 'Variance'])
    #plt.ylim(0.52, 0.62)
    
    
#%%
    
    sim = 100
    
    df = pd.DataFrame({'SIM':range(0,sim)})
    
    fig = plt.figure()
    
    for i,sm0 in enumerate(np.linspace(0,1.0,11)):
        # Run model 
        result = annual_sm(sim,sm0)
        # Extract SM values for first day
        doy = result.iloc[0,:].reset_index(drop=True)
        df['SM_' + str(i)] = doy
        
        # Plot
        plt.plot(df.index,doy)
    
    plt.xlabel('Simulation')
    plt.ylabel('Soil moisture')
    plt.legend(labels=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],bbox_to_anchor=(1.02,1), loc="upper left")
    
    


#%%
if __name__ == "__main__":
    main()
