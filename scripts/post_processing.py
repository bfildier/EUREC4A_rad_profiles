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

#Constants

epsilon = 0.6223 # Ratio of molar mas of water to dry air
p_0 = 1.e5 #Pa
kappa = 2/7 #Poisson constant
c_p = 1005 #J/kg/K
Ra = 286.9
g = 10 #m.s^-2
day_to_s = 86400 #s


def concat(input_dir):
    
    all_rad_path = glob.glob(os.path.join(input_dir,'*_rrtmgp.nc'))
    all_rad_path.sort()

    all_rad_files = []

    for rad_path in all_rad_path:

        sonde = xr.open_dataset(rad_path)
        print("")
        print(rad_path)
        
        if (compute_q_rad == True):
            sonde = calculateMixingRatio(sonde)
            sonde = calculateDensity(sonde)

           #Calculate altitude of level using hydrostatic approximation

            n_play = len(sonde["play"])
            n_plev = len(sonde["plev"])

            zlay = np.zeros(n_play)
            zlev = np.zeros(n_plev)

            zlev[n_plev-1] = 0

            for k in range(n_plev-1,0,-1):
                zlev[k-1] = zlev[k] +  1/(sonde.rho.values[k-1]*g)*(sonde.plev.values[k] - sonde.plev.values[k-1])

            zlay[n_play-1] = 0.5*(zlev[n_plev-1] + zlev[n_plev-2])

            for k in range(n_play-1,0,-1):
                rho_lev = 0.5*(sonde.rho.values[k-1]+sonde.rho.values[k])
                zlay[k-1] = zlay[k] + 1/(rho_lev*g)*(sonde.play.values[k] - sonde.play.values[k-1])

            sonde["zlev"] = (["plev"],zlev)
            sonde["zlay"] = (["play"], zlay)
                
            sonde = calculateQrad(sonde)
    
                
            sonde = sonde.swap_dims({'play':'zlay'}).swap_dims({'plev':'zlev'}).\
            interp(zlev=np.arange(0,10005,10)).interp(zlay=np.arange(5,10000,10))    


            
        else:
            
            sonde = sonde.interp(plev=np.arange(101500,50,-100)).interp(play=np.arange(101450,50,-100))
        
        all_rad_files.append(sonde)


    all_rad_files = xr.concat(all_rad_files, dim="launch_time")
    
    return all_rad_files


def calculateMixingRatio(all_sondes):    
    
    all_sondes = all_sondes.assign(mr=all_sondes["h2o"]*epsilon)   #volume to mass mixing ratio
    
    return all_sondes
    

def calculateDensity(sonde):
    
    #Formula rho=pres/RaT*(1+x)/(1+x*Rw/Ra)

    sonde = sonde.assign(rho=sonde["play"]/(Ra*sonde["tlay"])*(1+sonde["mr"])/
                                                           (1+1.609*sonde["mr"]))    
    return sonde

def calculateQrad(sonde):
    
    plev_size = len(sonde.plev)
    play_size = len(sonde.play)
    
    flux = np.zeros(plev_size)
    flux_lw = np.zeros(plev_size)
    flux_sw = np.zeros(plev_size)
    delta_zlev = np.zeros(plev_size)
    
    q_rad = np.zeros(play_size)
    q_rad_lw = np.zeros(play_size)
    q_rad_sw = np.zeros(play_size)
    
    for k in range(plev_size):
        flux_lw[k] = sonde.lw_up[k] - sonde.lw_dn[k]
        flux_sw[k] = sonde.sw_up[k] - sonde.sw_dn[k]
        flux[k] = flux_lw[k] + flux_sw[k]
        
    for k in range(plev_size-1):
        delta_zlev[k] = sonde.zlev[k+1] - sonde.zlev[k]
        
    for k in range(play_size):
        q_rad_lw[k] = - 1/(sonde.rho[k]*c_p)*day_to_s*(flux_lw[k+1]-flux_lw[k])/delta_zlev[k]
        q_rad_sw[k] = - 1/(sonde.rho[k]*c_p)*day_to_s*(flux_sw[k+1]-flux_sw[k])/delta_zlev[k]
        q_rad[k] = - 1/(sonde.rho[k]*c_p)*day_to_s*(flux[k+1]-flux[k])/delta_zlev[k]

    
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
    parser.add_argument("--comp_qrad", type=bool, default=False,
                        help="to compute the radiative cooling for the sondes")
    args = parser.parse_args()

    # Generalize this
    input_dir = args.in_dir
    output_dir = args.out_dir
    compute_q_rad = args.comp_qrad

all_rad_files = concat(input_dir)
all_rad_files.to_netcdf(output_dir+"all_rad_profiles_radiosondes.nc")