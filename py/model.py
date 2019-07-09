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
        from numpy import zeros
        self.soil = kwargs.pop('soil')
        self.crop = kwargs.pop('crop')
        self.climate = kwargs.pop('climate')
        self.n_days = len(self.climate.rainfall)

        # Set maximum soil water content (mm):
        self.nZr = self.soil.n * self.crop.Zr # [mm]  

        # Pre-allocate arrays
        self.R = self.climate.rainfall
        self.s = zeros(self.n_days)     # [0-1]
        self.ET = zeros(self.n_days)
        self.E = zeros(self.n_days)
        self.T = zeros(self.n_days)
        self.L = zeros(self.n_days)
        self.dsdt = zeros(self.n_days)  # [mm/day]
        
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
                    self.L[t]  = (self.s[t+1] - 1)*self.nZr # [mm/day]
                    self.dsdt[t] = self.dsdt[t] - self.L[t] # [mm/day]
                    self.s[t+1] = 1
                if self.s[t+1] < 0:
                    self.s[t+1] = 0
                
                # 4. Determine leakage loss:
                if self.s[t+1] > self.soil.sfc:
                    # Calculate leakage in units of mm/day
                    self.L[t] = self.L[t] + (self.s[t+1] - self.soil.sfc)*self.nZr # [mm/day]
                    # Update dsdt for leakage loss
                    self.dsdt[t] = self.dsdt[t] - (self.s[t+1] - self.soil.sfc)*self.nZr # [mm/day]
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
            'ET':self.E + self.T,
            'T':self.T,
            'L':self.L,
            'dsdt':self.dsdt,

        })
