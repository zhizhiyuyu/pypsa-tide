
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad


def read_area_elevation_curve(file, depth_correction  = 4, conversion=1.0):
    """
    Reading depth-area curve.
    :param file: input_file with array of format (h, area)
    :param depth_correction : amplitude correction in m (depending on datum)
    :return: h-area curve function
    """

    depth_area = np.load(file)
    return interp1d(depth_area[:,0] - depth_correction,depth_area[:,1] * conversion)


def lagoon_volume_calculation(file, depth_correction=4, elevation=0):

    depth_area = np.load(file)
    f = read_area_elevation_curve(file,depth_correction=depth_correction, conversion=1e6)  # Conversion - area is in km2 need to be consistent with units and convert to m2
    volume = quad(f,depth_area[0,0],elevation, limit=50,)
    return volume

# Let's assume a constant area over the a fixed depth of 10 m  (Answer should be 1e7)
a = lambda x: 1e6
volume = quad(a,-10,0, maxp1=40, limlst=40)

# Now let's calculate the volume from the area swansea curve
vol = lagoon_volume_calculation('../inputs/area_swansea.npy')

# Printing the volume calculated using scipy for the idealised case and based on the depth area curve
print(volume, vol, volume[0]-vol[0])
