import numpy as np
from netCDF4 import Dataset, num2date # to work with NetCDF files
from os.path import expanduser
import matplotlib.pyplot as plt
home = expanduser("~") # Get users home directory
from scipy import stats
plt.rcParams.update({'font.size': 35})

import xarray as xr
import glob, os
import numpy as np
import matplotlib.pyplot as plt      
import datetime
import math
import metpy
from metpy import calc as mpcalc
import pandas as pd

Rd=287
g=9.8

direc = "/media/ludo/DATA/google-drive/Thèse/EUREC4a/github/Input/ERA"
q_path = os.path.join(direc, "q_2020_01_ERA5_Barbados.nc")
t_path = os.path.join(direc, "ta_2020_01_ERA5_Barbados.nc")


#We select the profiles corresponding to the time and location of the sounding
q = xr.open_dataset(q_path).isel(time=1, latitude=1, longitude=1)
t = xr.open_dataset(t_path).isel(time=1, latitude=1, longitude=1)


one_profile = q
one_profile["ta"] = t["ta"]

#We need the following data, corresponding to the highest level measured by the dropsonde
alt_max_sonde = 10000 #in m
pres_max_sonde = 250 #in hPa
temp_sonde = 225 #in K (one could use directly the temperature instead)

#In ERA profile, we keep only pressure values inferior to the highest pressure measured by the sounding
one_profile_above = one_profile.where(one_profile.level<pres_max_sonde, drop=True)
length = len(one_profile_above.level)


#array where we will store the altitude
altitude_above = []

#We initialize this array
first_level_above = one_profile.isel(level=length-1)
first_alt_above = alt_max_sonde + Rd/g*((first_level_above.ta.values + temp_sonde)/2\
                                        *math.log(pres_max_sonde/first_level_above.level.values))
altitude_above.insert(0, first_alt_above)

#We calculate the altitude of each level using hydrostatic equation
for i in range(length-1, 1, -1):
    profile_1 = one_profile.isel(level=i)
    profile_2 = one_profile.isel(level=i-1)
    
    new_alt = altitude_above[0] + Rd/g*((profile_1.ta.values + profile_2.ta.values)/2
                                        *math.log(profile_1.level.values/profile_2.level.values))
    
    altitude_above.insert(0, new_alt)
    
print(altitude_above)

