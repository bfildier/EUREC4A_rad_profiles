#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 12:34:19 2020

@author: bfildier
"""


import xarray as xr
import numpy as np
import os
from scipy.interpolate import interp2d


input_dir = '/Users/bfildier/Data/EUREC4A/ERA'
q_example = xr.open_dataset(os.path.join(input_dir,'q_2020_01_02_ERA5_hourly_Barbados.nc'))

def interpScalar(values,lat_coord,lon_coord,lat,lon):
    
    # We can think of changing option 'kind=' to something different than a linear interpolation
    return interp2d(lat_coord,lon_coord,values)(lat,lon)

def interpProfile(values,lat_coord,lon_coord,lat,lon):

    Nz = values.shape[0]
    profile_out = np.nan*np.zeros((Nz,))
    
    for i_z in range(Nz):
        profile_out[i_z] = interpScalar(values[i_z],lat_coord,lon_coord,lat,lon)
    
    return profile_out

lon_coord = q_example.longitude.values
lat_coord = q_example.latitude.values

# TEST

print('-- test scalar interpolation')
print(' interpolated value:')
print(interpScalar(q_example.q[0].values[0],lat_coord,lon_coord,16.9,300.1))
print(' neighboring values:')
print(q_example.q[0].values[0,:2,:2])
print()
print('-- test profile interpolation')    
print(interpProfile(q_example.q[0].values,lat_coord,lon_coord,16.9,300.1))