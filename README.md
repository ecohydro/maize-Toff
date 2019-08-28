# maize-Toff
Repo for modeling analysis of maize yield variability and tradeoffs between yield and crop failure.

## Research questions
1. What are the tradeoffs between yield (i.e. cultivar or variety selection), planting date, and crop failure?
2. To what extent does temperature and precipitation variability drive variability in crop production?

## Set-up
1. ```cd maize-Toff```
2. ```conda env create -f environment.yml -n maize-Toff```
3. ```conda activate maize```
4. ```jupyter notebook```

Note: To update dependecies in an existing environment, use `conda env update --file environment.yml` after step four.

## Directory structure

### farm/
* where models are stored

### data/
* contains CETRAD rainfall data, maize variety info, survey info...

### output/
* exported figures, results, etc.

Draft of manuscript [here](https://drive.google.com/drive/folders/1uRBOhii9ca5uGy2ntIhcAtsCL4PKBGq-?usp=sharing).

## How to use this model

### Generate a Soil object

The soil object contains all the necessary parameters to specify soil hydraulic properties and texture. The easiest
way to generate a soil instance is to use one of the standard USDA soil texture types: `Sand`, `Loamy Sand`, `Sandy Loam`, `Silt Loam`, `Loam`, `Sandy Clay Loam`, 
`Clay Loam`, `Sandy Clay`, `Silty Clay`, `Clay`, and `Sandy Silty Loamy Clay`. Just kidding. That last one is totally not a soil type.

```python
# Creates a sand soil object:
soil = Soil('sand') 
```
You can also specify your own parameters:
```python

param_dict = {
    'b': 11.4,
    'Psi_S_cm': 40.5,   # saturated water tension, cm
    'Psi_l_cm': 24.3,   # leakage water tension, cm
    'n': 0.482,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
    'Ks': 0.0077,   # saturated hydraulic conductivity, cm/min
    'S': 0.268      # sorptivity, cm/min^1/2  
}

custom_soil = Soil(params=param_dict)
```

Of the parameters, only `b`, `Psi_S_cm`, `Psi_l_cm`, `n`, and `Ks` are required. The model doesn't use `S` right now.

### Generate a climate object

The climate object has information on maximum evapotranspiration, as well as a time series of rainfall.
The rainfall timeseries is generated stochastically using
the parameters of storm depth, frequency, and the length of the season. The storm depth (`alpha_r`) specifies the average daily storm depth, assuming that daily storm depths are drawn from an exponential distribution. The frequency (`lambda_r`) is best described as the daily probability of rainfall. The length of the season ends up setting the timescale of the simulation in days.

```python
# Average storm depth [mm]
alpha_r=10.0    

# Frequency of storms [day^-1]
lambda_r=0.3   

 # Length of rainy season [days]
t_seas=180     

# Maximum evapotranspiration [mm/day]
ET_max=6.5      

climate = Climate(
    alpha_r=alpha_r,
    lambda_r=lambda_r,
    t_seas=t_seas,
    ET_max=ET_max)

```

### Generate a plant

The plant object requires the most parameters to generate, and some of these parameters depend on the soil in which the plant is growing. In addition, the plant object can be subclassed into specific plant types to allow for varying structures and function. In this simulation, we are using the `Crop` subclass, which is initialized with the minimum following parameters:

* `kc_max` The maximum crop coefficient [`dimensionless`], which is a scale factor applied to `ET_max` to determine the maximum rate of plant transpiration, `T_max`.

* `LAI_max` The maximum crop leaf area [`m^2/m^2`].

* `T_max` The maximum rate of crop water use in [`mm/day`]

* `soil` A soil object that specifies the soil that this crop is growing in.

```python

crop = Crop(kc_max=1.2, LAI_max=2.0, T_max=4.0, soil=soil)

```

Additional parameters that _should_ be specified are:

* `Zr` Plant rooting depth [`mm`]. Default is 500mm.
* `sw_MPa` Plant wilting point [`MPa`]. Default is -1.5MPa
* `s_star_MPa` Water potential of maximum water use [`MPa`]. Default is -0.2 MPa.

If these parameters are not provided, default values are inherited from the `Plant` class.

### Initialize and run the model

Combining the prior three steps, we can create a model instance and run the model:

```python

# Import the necessary objects:
from farm.climate import Climate
from farm.soil import Soil
from farm.plant import Crop
from farm.model import CropModel

# Make the things
climate = Climate() # uses default climate values
soil = Soil('sand')
crop = Crop(kc_max=1.2, LAI_max=2.0, T_max=4.0, soil=soil)
model = Model(crop=crop,soil=soil,climate=climate)

# RUN IT.
model.run() # TADA!

```
### Getting model output 

The `model.output()` function returns all the simulation output structured as a single pandas DataFrame. The frame has the following columns:

* `kc` Time series of daily crop coefficients. 
* `LAI` Time series of crop LAI
* `R` Time series of daily rainfall [mm]
* `s` Time series of daily relative soil moisture [0-1]
* `E` Time series of daily soil evaporation [mm]
* `T` Time series of daily plant transpiration [mm]
* `L` Time series of soil leakage loss [mm]
* `Q` Time series of surface runoff [mm]
* `dsdt` Time series of changing relative soil moisture

To plot any of this data, simply use the `.plot()` command:

```python

output = model.output()

# Plots a time series of simulated evapotranspiration:
output['ET'].plot()

# Plots a time series of simulated relative soil moisture
output['s'].plot()

```
