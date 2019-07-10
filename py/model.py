"""
Name:           model.py
Compatibility:  Python 3.7
Description:    Ecohydro model to track evaporation and transpiration

URL:            https://github.com/ecohydro/maize-Toff

Requires:       pandas; matlotlib; numpy; math

AUTHOR:         Natasha Krell & Kelly Caylor
ORGANIZATION:   U.C. Santa Barbara
Contact:        nkrell@ucsb.edu

"""
#%% Define the CropModel Class
from pandas import DataFrame
from numpy import zeros
        
class CropModel():
    """ Defines an ecohydrological model class for use with crops.

    Usage:

        model = CropModel(
            soil=soil,
            climate=climate,
            crop=crop
        )

    """
    def __init__(self, *args, **kwargs):
        self.soil = kwargs.pop('soil')
        self.crop = kwargs.pop('crop')
        self.climate = kwargs.pop('climate')
        self.n_days = len(self.climate.rainfall)

        try:
            # Set the nZr using the soil's function.
            self.nZr = soil.set_nZr(crop)
        except:
            # Set maximum soil water content (mm):
            self.nZr = self.soil.n * self.crop.Zr # [mm]  

        # Pre-allocate arrays
        self.R = self.climate.rainfall  # [mm/day]
        self.s = zeros(self.n_days)     # [0-1]
        self.ET = zeros(self.n_days)    # [mm/day]
        self.E = zeros(self.n_days)     # [mm/day]
        self.T = zeros(self.n_days)     # [mm/day]
        self.L = zeros(self.n_days)     # [mm/day]
        self.Q = zeros(self.n_days)     # [mm/day]

        # Note: We calculate dsdt in mm/day to make it easier 
        # to handle water balance wrt to other terms.
        self.dsdt = zeros(self.n_days)  # [mm/day]
        
        self.LAI = zeros(self.n_days)
        self.ET_max = zeros(self.n_days)
        self.T_max = zeros(self.n_days)
        self.kc = zeros(self.n_days)

        # Set initial conditions:
        self.s[0] = 0.3     # relative soil moisture, [0-1]
        _s = self.s[0]      # intermediate soil moisture used during
                            # model time step calculations.
        _dsdt = 0           # intermediate soil moisture change used
                            #  during model time step calculations.
        
    def run(self):
        for t in range(self.n_days):
            try:
                # 0. Update the crop coefficient
                # TODO: Edit Crop class to make this dynamic.
                self.kc[t] = self.crop.calc_kc(t)
                self.LAI[t] = self.crop.calc_LAI(t)

                # 2. Calculate ET terms
                self.T[t] = self.crop.calc_T(self.s[t],t)   # mm/day
                self.E[t] = self.climate.calc_E(
                    self.s[t],
                    t,
                    plant=self.crop,
                    soil=self.soil) # mm/day
                self.ET[t] = self.T[t] + self.E[t]
                
                # 1. Update Soil Moisture Water Balance (Part 1)

                """ Note:

                The model does water balance in two parts. 
                The first part determines the effect of any Rainfall 
                and Evapotranspiration on the relative soil moisture.
                The input of rainfall could easily cause the soil moisture
                to get above either the field capacity (soil.sfc) or saturation
                (i.e. s > 1). In either case, two additional fluxes would kick in.

                Therefore, we first see how Rainfall and ET affect s, and update a
                temporary s value, _s, with the temporary dsdt value, _dsdt.

                _dsdt = R(t) - E(s,t)
                _s = s(t) + _dsdt

                Now we can use the value of _s to see if leakage or runoff
                occur in this time step (Step 3 and Step 4). 
                
                Then, once we have values of L(t) and Q(t) based on _s,
                we update the dsdt and s values using our master equations:

                dsdt(t) = R(t) - E(s,t) - L(s,t) - Q(s,t)
                s(t+1) = s(t) + dsdt

                All this _s and _dsdt stuff is an abstraction that allows us to 
                ensure that we handle water balance in a way that makes sense.

                Deal with it.

                """
                # Create temporary s and dsdt value for 
                # use within timestep calculations:
                _dsdt = self.R[t] - self.ET[t]         # mm/day
                _s = self.s[t] + _dsdt/self.nZr
                
                # 3. Determine saturation excess flow
                # Note: Use the temporary (intermediate) s value 
                # for this calculation rather than s[t]
                self.Q[t] = self.soil.calc_Q(
                    _s, units='mm/day')
                
                # 4. Determine leakage loss:
                # Note: Use the temporary (intermediate) s value 
                # for this calculation rather than s[t]
                self.L[t] = self.soil.calc_L(
                    _s, units='mm/day')

                # 5. Update Soil Moisture Water Balance (Part 2)
                self.dsdt[t] = self.R[t] - self.ET[t] - self.Q[t] - self.L[t]
                self.s[t+1] = self.s[t] + self.dsdt[t]/self.nZr
                # print(
                #     f"Time: {t}\t s(t):{self.s[t]:.3f}\t"
                #     f"dsdt:{self.dsdt[t]:.3f}\t s(t+1):{self.s[t+1]:.3f}"
                #     f"Q[t]:{self.Q[t]:.3f}\t L[t]:{self.L[t]:.3f}"
                # )
            except IndexError:
                print("DONE. At end of simulation, timestep {t}".format(t=t))
                break
    
    def output(self):
        return DataFrame({ 
            'kc':self.kc,
            'LAI':self.LAI,
            'R':self.R,
            's':self.s,
            'E':self.E,
            'ET':self.E + self.T,
            'T':self.T,
            'L':self.L,
            'dsdt':self.dsdt,

        })
