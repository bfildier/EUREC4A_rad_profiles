import numpy as np
from netCDF4 import Dataset, num2date # to work with NetCDF files
import xarray as xr
import glob, os
import numpy as np
import matplotlib.pyplot as plt      
import datetime
import metpy
import argparse
import pytz
from metpy import calc as mpcalc

#Constants

epsilon = 0.6223 # Ratio of molar mas of water to dry air
p_0 = 1.e5 #Pa
kappa = 2/7 #Poisson constant
c_p = 1005 #J/kg/K
Ra = 286.9
g = 9.8 #m.s^-2
day_to_s = 86400 #s

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def concat(input_dir):
    
    all_rad_path = glob.glob(os.path.join(input_dir,'rrtmgp_*.nc'))
    all_rad_path.sort()

    all_rad_files = []

    for rad_path in all_rad_path:

        sonde = xr.open_dataset(rad_path)
        print("")
        print(rad_path)
        
        if (compute_q_rad == True):
            sonde = calculateMixingRatio(sonde)
            sonde = calculateDensity(sonde)
            
            sonde = calculateQrad(sonde)
                
            if args.interp:
                sonde = sonde.swap_dims({'play':'zlay'}).swap_dims({'plev':'zlev'})
                sonde = sonde.interp(zlev=np.arange(0,10005,10)).interp(zlay=np.arange(5,10000,10)) 
        
        else:
            
            if args.interp:
                sonde = sonde.interp(plev=np.arange(101500,50,-100)).interp(play=np.arange(101450,50,-100))
        
        all_rad_files.append(sonde)


    all_rad_files = xr.concat(all_rad_files, dim="launch_time")
    
    return all_rad_files


def calculateMixingRatio(all_sondes):    
    
    all_sondes = all_sondes.assign(mr=all_sondes["h2o"]*epsilon)   #volume to mass mixing ratio
    
    return all_sondes
    

def calculateDensity(sonde):
    
    #Formula rho=pres/RaT*(1+x)/(1+x*Rw/Ra)
    
    sonde["play"].attrs['units'] = 'Pa'
    sonde["tlay"].attrs['units'] = 'K'

    rho = mpcalc.density(sonde["play"], sonde["tlay"], sonde["mr"])
    
    sonde["rho"] = (["play"], rho.magnitude)
                  
    return sonde

def calculateQrad(sonde):
        
    flux_lw = sonde.lw_up[:] - sonde.lw_dn[:]
    flux_sw = sonde.sw_up[:] - sonde.sw_dn[:]
    flux = flux_lw[:] + flux_sw[:]  
      
    play_size = len(sonde.play)
    delta_plev = np.zeros(play_size)
    for k in range(play_size):
        delta_plev[k] = sonde.plev.values[k+1] - sonde.plev.values[k]

    q_rad = np.zeros(play_size)
    q_rad_lw = np.zeros(play_size)
    q_rad_sw = np.zeros(play_size)
        
    for k in range(play_size):  
        q_rad_lw[k] = 1/(c_p)*g*day_to_s*(flux_lw[k+1]-flux_lw[k])/delta_plev[k]
        q_rad_sw[k] = 1/(c_p)*g*day_to_s*(flux_sw[k+1]-flux_sw[k])/delta_plev[k]
        q_rad[k] = 1/(c_p)*g*day_to_s*(flux[k+1]-flux[k])/delta_plev[k]

    
    sonde["q_rad"] = (("play"), q_rad)
    sonde["q_rad_lw"] = (("play"), q_rad_lw)
    sonde["q_rad_sw"] = (("play"), q_rad_sw)
    
    return sonde



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Concatenate all radiative profiles into one single file")
    parser.add_argument("--in_dir", type=str, default="../output/rad_profiles",
                        help="Directory where the profiles to concatenate are saved")
    parser.add_argument("--out_dir", type=str, default="../output/",
                        help="Directory where the output file should be saved")
    parser.add_argument("--comp_qrad", type=str2bool, nargs='?',
                        const=True, default=False,
                        help="to compute the radiative cooling for the sondes")
    parser.add_argument("--interp",type=str2bool, nargs='?',
                        const=True, default=True,
                        help="interpolate on finer grid (default = True)")
    args = parser.parse_args()

    # Generalize this
    input_dir = args.in_dir
    output_dir = args.out_dir
    compute_q_rad = args.comp_qrad

    print('comp_qrad =',args.comp_qrad)
    print('interp =',args.interp)

    all_rad_files = concat(input_dir)
    all_rad_files.to_netcdf(os.path.join(output_dir,"rad_profiles.nc"), mode="w")