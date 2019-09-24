# import objects
from farm import Climate
from farm import Soil
from farm import Crop
from farm import CropModel

# initialize objects
climate = Climate() # uses default climate values
soil = Soil('sand')
crop = Crop(soil=soil) 
soil.set_nZr(crop)
model = CropModel(crop=crop,soil=soil,climate=climate)

model.run()

# write tests
def test_sh():
	assert soil.sh <= crop.sw, "Should be true"

def test_star():
	assert crop.s_star <= soil.sfc, "Error: s_star is too high: it should be less than field capacity"


if __name__ == "__main__":
    test_sh()
    test_star()
    print("Everything passed")