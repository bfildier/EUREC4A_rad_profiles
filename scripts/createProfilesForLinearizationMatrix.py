#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 16:17:28 2020

@author: bfildier
"""

# load modules

import matplotlib.pyplot as plt
import argparse
import os,glob
from datetime import datetime

import numpy as np
import xarray as xr



#%% Define directories

# current dir
wdir = os.path.dirname(os.path.realpath(__file__))
# code dir
cdir = os.path.join(os.path.dirname(wdir),'code')
# output dir
odir = os.path.join(os.path.dirname(wdir),'output/linear_model')
odir_rad_profiles = os.path.join(odir,'rad_profiles')
os.makedirs(odir,exist_ok=True)
os.makedirs(odir_rad_profiles,exist_ok=True)

#%% Read arguments

parser = argparse.ArgumentParser(description="Calculate radiative profiles corresponding to individual T and q perturbations at each level to be used to linearize the problem")
parser.add_argument("--sonde_file", type=str, default="/Users/bfildier/Data/EUREC4A/merged/sondes/all_sondes.nc",
                    help="Name of sonde file")
args = parser.parse_args()

#%% Do calculation

#%%    # 0. calculate mean/std of sonde data and save

all_sondes = xr.open_dataset(args.sonde_file)
# Create latitude/longitude variables 
lat_mean = all_sondes.latitude.mean(dim='launch_time')
lon_mean = all_sondes.longitude.mean(dim='launch_time')
# Calculate mean sounding
sondes_mean = all_sondes.mean(dim='launch_time')
sondes_std = all_sondes.std(dim='launch_time')

def addDimensionsAndSave(sondes,savename):

    # add launch_time dimension
    sondes = sondes.assign_coords(coords={'launch_time':datetime(year=2020,month=2,day=1)})
    sondes = sondes.expand_dims('launch_time')
    # add spatial coordinates and platform variables
    sondes = sondes.assign({'latitude':lat_mean,
                            'longitude':lon_mean,
                            'Platform':'all'})
    # save
    sondes.to_netcdf(os.path.join(odir,savename))

# save
addDimensionsAndSave(sondes_mean,'sondes_mean.nc')
addDimensionsAndSave(sondes_std,'sondes_std.nc') # for fun (or just in case)

#%%    # 1. combine mean profile with reference background

print(' ---- combine mean profile with reference profile')

combine_command = "cd %s; python combine_reference_and_sonde_profiles.py"%wdir+\
    " --sonde_file=%s --out_dir=%s"%(os.path.join(odir,'sondes_mean.nc'),
                                     odir)

os.system(combine_command)

#%%    # 2. Interpolate values onto a smaller grid (100 regular levels in pressure)

# reload combined mean profiles
sondes_mean_combined = xr.open_dataset(os.path.join(odir,'rrtmgp_0000_all.nc'))
# # temporary
# all_sondes_combined = xr.open_dataset(os.path.join(odir,'../rad_profiles_all_sondes.nc'))

# create new vertical grid
Np = 100 # number of levels
plev_new = np.linspace(10,1000,Np)*100 # in Pa
play_new = np.convolve(plev_new,[0.5,0.5],mode='valid')

# interpolate dataset onto new grid
var_dict = {}
for varid in list(sondes_mean_combined.data_vars):
    
    if varid in ['launch_time','platform','sfc_emis','sfc_alb','sfc_t','cos_sza']:
        
        # reassign directly
        var_dict[varid] = ([],sondes_mean_combined[varid])
        
    elif varid in ['tlay','h2o','zlay','co2','ch4','n2o','o3','o2','n2','co']:
        
        # interpolate onto play_new
        var_dict[varid] = (["play"],sondes_mean_combined[varid].interp(play=play_new))
        
    else:
        
        # interpolate onto plev_new
        var_dict[varid] = (["plev"],sondes_mean_combined[varid].interp(plev=plev_new))

# add coordinates
var_dict['play'] = (["play"],play_new)
var_dict['plev'] = (["plev"],plev_new)

# create new dataset with interpolated variables
sondes_mean_interp = xr.Dataset(var_dict)

#%%    # 3. generate separate files with pertubed T or q at individual levels
#        use sondes_std to define perturbation

# list of all sondes to feed the radiative transfer code
sondes_perturbed = [sondes_mean_interp]

# set perturbations magnitudes
perturb = {'tlay' : 1, # K
           'h2o': 0.001} # kg/kg

# apply
for varid in perturb.keys():
    
    for i_p in range(1,Np):
        
        # copy mean
        new_sonde = sondes_mean_interp.copy()
        # vary layer i_p by chosen perturbation
        new_values = new_sonde[varid].values.copy()
        new_values[-i_p] += perturb[varid] # replace item #i_p starting from end
        # reassign
        new_sonde[varid].values = new_values
        # store
        sondes_perturbed.append(new_sonde)
        
# save all profiles
for i,sonde in zip(range(len(sondes_perturbed)),sondes_perturbed):
    
    name = 'rrtmgp_%s.nc'%('{:>04}'.format(i))
    sonde.to_netcdf(os.path.join(odir_rad_profiles,name))
    
#%%    # 4. run rrtmgp on mean profile and on perturbed profiles

for i in range(2*Np+1):

    name = 'rrtmgp_%s.nc'%('{:>04}'.format(i))
    print(name)

    compute_rad_command = "cd %s; ./sonde_radiation %s"%(cdir,
                                                   os.path.join(odir_rad_profiles,name))

    os.system(compute_rad_command)

#%%    # 5. postprocessing
    
postprocessing_command = "cd %s; python post_processing.py --in_dir=%s --out_dir=%s --comp_qrad=True --interp=False"%\
    (wdir,odir_rad_profiles,odir)
    
os.system(postprocessing_command)


    