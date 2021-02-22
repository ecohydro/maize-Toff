"""
Name:           model.py
Compatibility:  Python 3.7
Description:    Ecohydro model to track evaporation and transpiration

URL:            https://github.com/ecohydro/maize-Toff

Requires:       pandas; matlotlib; numpy; math; logging

AUTHOR:         Natasha Krell & Kelly Caylor
ORGANIZATION:   U.C. Santa Barbara
Contact:        nkrell@ucsb.edu

"""
#%% Define the CropModel Class
from pandas import DataFrame
from numpy import zeros
import numpy as np
from .climate import Climate
import functools

import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
        
class CropModel():
    """ Defines an ecohydrological model class for use with crops.

    Usage:

        model = CropModel(
            soil=soil,
            climate=climate,
            crop=crop,
            planting_date=100
        )

    Default values:
        planting_date = 100 # Date Planted [Julian Day]

    """
    def __init__(self, soil=None, climate=None, crop=None, *args):
        """
        Initializes a crop model object.

        The init function requires the user to input crop, soil, and climate objects.
        It has a default planting date of 100, which is early April; a typical 
        value for the long rains cropping season in Kenya (Ray et al. 2015)

        """
        self.soil = soil
        self.crop = crop
        self.climate = climate
        #self.n_days = len(self.climate.rainfall)
        
        try:
            # Set the nZr using the soil's function.
            self.nZr = soil.set_nZr(crop)
        except:
            # Set maximum soil water content (mm):
            self.nZr = self.soil.n * self.crop.Zr # [mm]  


    def pre_allocate(self, n_days=None):

        # Pre-allocate arrays
        self.R = zeros(self.n_days)     # Will add in rainfall later
        self.s = zeros(self.n_days)     # [0-1]
        self.ET = zeros(self.n_days)    # [mm/day]
        self.I = zeros(self.n_days)     # [mm/day]
        self.E = zeros(self.n_days)     # [mm/day]
        self.T = zeros(self.n_days)     # [mm/day]
        self.L = zeros(self.n_days)     # [mm/day]
        self.Q = zeros(self.n_days)     # [mm/day]
        self.dos = zeros(self.n_days)   # integer

        # Note: We calculate dsdt in mm/day to make it easier 
        # to handle water balance with regards to other terms.
        self.dsdt = zeros(self.n_days)  # [mm/day]
        
        self.LAI = zeros(self.n_days)
        self.ET_max = zeros(self.n_days)
        self.T_max = zeros(self.n_days)
        self.kc = zeros(self.n_days)
        self.stress = zeros(self.n_days)        

    def run(self,
            do_output=False,
            s0=0.3,
            planting_date=100,
            t_before=60,
            t_after=7):
        """ Runs the ecohydro crop model.

            Inputs:
                planting_date: This is the julian day of 
                    year that the crop is planted.

                s0: This is the initial soil moisture value 
                    for the simulation

                t_before: This is the time in days before planting date
                    that we start the model.

                t_after: This is the duration after harvest date that we
                    run the model.

            We initialize the model run to start on DOY planting_date - t_before.
            The model then runs for a total of t_before + crop.lgp + t_after days

        """
        self.planting_date = planting_date
        self.n_days = t_before + self.crop.lgp + t_after
        self.pre_allocate()

        self.doy_start = self.planting_date - t_before
        # Make sure we don't back into the prior year.
        if self.doy_start <= 0:
            self.doy_start = 365 + self.doy_start

        self.dos_end = self.planting_date + self.crop.lgp + t_after

        # Force doy to be in [1,365]:
        doy = np.arange(self.doy_start, self.doy_start + self.n_days)
        while (doy - 365 > 0).any() == True:
            doy = doy - 365 * ((doy - 365) > 0)

        self.doy_end = doy[-1:] # last doy of year of simulation.

        for t in range(self.n_days):
            self.R[t] = self.climate.rainfall[doy[t]-1]

        self.doy = doy

        self.s[0] = s0     # relative soil moisture, [0-1]
        _s = self.s[0]      # intermediate soil moisture used during
                            # model time step calculations.
        _dsdt = 0           # intermediate soil moisture change used
                            # during model time step calculations.
        
        dos = 0   # initialize day of season to zero.
        planted = False # flag to determine if the season has started.

        for t in range(self.n_days):
            try:        
                # Determine the day of season.
                if self.doy[t] == self.planting_date: 
                    planted = True
                if planted == True:
                    dos = dos + 1 # t_seas init as zero
                # 0. Update the crop coefficient and day of season
                self.dos[t] = dos
                self.kc[t] = self.crop.calc_kc(dos)
                self.LAI[t] = self.crop.calc_LAI(self.kc[dos])
                self.stress[t] = self.crop.calc_stress(self.s[t])

                # 1. Calculate Q
                self.Q[t] = self.soil.calc_Q(self.s[t], units='mm/day')


                # 2. Update Soil Moisture Water Balance (Part 1)

                """ Note:

                The model does water balance in two parts. 
                The first part determines the effect of any Rainfall 
                and Evapotranspiration on the relative soil moisture.
                The input of rainfall could easily cause the soil moisture
                to get above either the field capacity (soil.sfc) or saturation
                (i.e. s > 1). In either case, two additional fluxes would kick in.

                Therefore, we first see how Rainfall and runoff affect s, and update 
                a temporary s value, _s, with the temporary dsdt value, _dsdt.

                _dsdt = R(t) - Q(s)
                _s = s(t) + _dsdt

                Now we can use the value of _s to determine leakage and ET
                in this time step (Step 3 and Step 4). Note: we allow both leakage
                and ET to start from the same _s value, which overestimates ET,
                rather than underestimating ET (as in a bucket model).
                
                Then, once we have values of L(t) and ET(t) based on _s,
                we update the dsdt and s values using our master equations:

                dsdt(t) = R(t) - Q(s,t) - E(s,t) - L(s,t)
                s(t+1) = s(t) + dsdt

                All this _s and _dsdt stuff is an abstraction that allows us to 
                ensure that we handle water balance in a way that makes sense.

                Deal with it.

                """
                # Create temporary s and dsdt value for 
                # use within timestep calculations:
                rainfall = max(self.R[t] - self.crop.calc_I(),0)
                _dsdt = rainfall - self.Q[t]         # mm/day
                _s = self.s[t] + _dsdt/self.nZr
                
                # 3. Calculate ET terms.
                # Note: Use the temporary (intermediate) s value
                # for this calculation rather than s[t].
                self.I[t] = min(self.crop.calc_I(),self.R[t])
                self.T[t] = self.crop.calc_T(
                    _s, LAI=self.LAI[t])   # mm/day
                self.E[t] = self.climate.calc_E(
                    _s, LAI = self.LAI[t], 
                    sh = self.soil.sh)
                self.ET[t] = self.T[t] + self.E[t]

                # 4. Determine leakage loss:
                # Note: Use the temporary (intermediate) s value 
                # for this calculation rather than s[t]
                self.L[t] = self.soil.calc_L(
                    _s, units='mm/day')

                # 5. Update Soil Moisture Water Balance (Part 2)
                self.dsdt[t] = rainfall - self.Q[t] - self.ET[t] - self.L[t]
                self.s[t+1] = self.s[t] + self.dsdt[t]/self.nZr
                # print(
                #     f"Time: {t}\t s(t):{self.s[t]:.3f}\t"
                #     f"dsdt:{self.dsdt[t]:.3f}\t s(t+1):{self.s[t+1]:.3f}"
                #     f"Q[t]:{self.Q[t]:.3f}\t L[t]:{self.L[t]:.3f}"
                # )
            except IndexError:
                #print(f"DONE. At end of simulation, timestep {t}")
                logger.info('logging is easier than I was expecting')
                break

        if do_output:
            return self.output()

    def output(self):
        return DataFrame({ 
            'kc':self.kc,
            'LAI':self.LAI,
            'stress':self.stress,
            'R':self.R,
            's':self.s,
            'I':self.I,
            'E':self.E,
            'ET':self.E + self.T,
            'T':self.T,
            'L':self.L,
            'dsdt':self.dsdt,
            'dos':self.dos,
            'doy':self.doy
        })