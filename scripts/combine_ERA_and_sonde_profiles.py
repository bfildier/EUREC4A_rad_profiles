import pysolar.solar as pysolar
import datetime
import argparse
import pytz
import xarray as xr
import numpy as np
import os
import pytz
import netCDF4
from metpy import calc as mpcalc
from scipy.interpolate import interp2d
import math


# Define parameters

hPaToPa = 100.                                          # Convert from hPa to Pa
PatohPa = 1 / 100                                       # Convert from Pa to hPa
epsilon = 0.6223                                        # Ratio of molar mass of water to dry air
CtoK    = 273.15                                        # Celsius to Kelvin
gtokg   = 1.e-3                                         # g/kg to kg/kg
ghgs    = ["co2", "ch4", "n2o", "o3", "o2", "n2", "co"] # greenhouse gases
Rd=287
g=9.8

def combine_sonde_and_background(all_sondes_file, background_file, ERA_dir, deltaP=100, sfc_emis=.98, sfc_alb=0.07, mu0=1., ghgs=ghgs):

    # Background sounding
    back_tropical_atm = xr.open_dataset(background_file)
    
    # add filenames to Argparse
    ta_ERA5 = xr.open_dataset(os.path.join(ERA_dir, 'ta_2020_01_02_ERA5_hourly_Barbados.nc'))
    SST_ERA5 = xr.open_dataset(os.path.join(ERA_dir, 'SST_2020_01_02_ERA5_hourly_Barbados.nc'))
    q_ERA5 = xr.open_dataset(os.path.join(ERA_dir, 'q_2020_01_02_ERA5_hourly_Barbados.nc'))
    
    SST_ERA5 = SST_ERA5.bfill(dim="longitude")
    
    all_sondes = xr.open_dataset(all_sondes_file)
    all_profiles = []
    
    number_sondes = len(all_sondes.launch_time)

    for i in range(number_sondes):
    # for i in range(100,106):

        alt_var = "height"
        p_var = 'pressure'
        q_var = 'specific_humidity'
        t_var = 'temperature'
        #rh_var = 'relative_humidity'
        
        sonde = all_sondes.isel(launch_time = i).dropna(dim=alt_var,\
                               subset=[alt_var, p_var,q_var,t_var],\
                                                   how="any")
        
        # minimum altitude to which radiosonde ascends, or from which dropsonde is released
        #maybe add these three arguments to arg parse
        size = 10
        min_alt = 200
        max_alt = 3000
        
        if (sonde[alt_var].values.size < size or sonde[alt_var].min() > min_alt or sonde[alt_var].max() < max_alt ):
            print("The sonde is empty ")
            sonde.close()            
            
        else:
            
            #select the corresponding ERA file 
            
            #select lat, lon and date of the sonde
            lat = sonde.latitude.dropna(dim=alt_var).values[0]
            lon = sonde.longitude.dropna(dim=alt_var).values[0]
            if (lat > 20 or lat < 4 or lon < -65 or lon > -50):
                print("outside SST bounds")
                sonde.close()
            else:
            
                fmt = '%Y-%m-%dT%H:%M:%S'
                date = datetime.datetime.strptime(str(sonde.launch_time.values).split('.')[0], fmt)
                #we add the timezone not for now, but for later when we will calculate the solar zenith angle
                date_UTC = pytz.timezone('UTC').localize(date)

                #Select closest SST, q and t profiles
                q_ERA5_sel = q_ERA5.sel(time=date_UTC, method="nearest")
                ta_ERA5_sel = ta_ERA5.sel(time=date_UTC, method="nearest")
                SST_ERA5_sel = SST_ERA5.sel(time=date_UTC, method="nearest")
                level = q_ERA5_sel.level.values*hPaToPa       #we want the levels to be in Pa                                         
                               
                #Select lat lon coordinates of the new grid
                #Maybe the same for q, ta et SST? In this case, keep only one of the following
                lon_coord_q = q_ERA5.longitude.values
                lat_coord_q = q_ERA5.latitude.values

                lon_coord_ta = ta_ERA5.longitude.values
                lat_coord_ta = ta_ERA5.latitude.values

                lon_coord_SST = SST_ERA5.longitude.values
                lat_coord_SST = SST_ERA5.latitude.values

                #Perform the interpolation

                #Note in longitude, we add 360 to convert to the correct form

                q_ERA5_interp = interpProfile(q_ERA5_sel.q.values,lat_coord_q,lon_coord_q,lat,lon+360)
                ta_ERA5_interp = interpProfile(ta_ERA5_sel.ta.values,lat_coord_ta,lon_coord_ta,lat,lon+360)
                SST_ERA5_interp = interpScalar(SST_ERA5_sel.sstk.values,lat_coord_SST,lon_coord_SST,lat,lon+360)
                

                #convert specific humidity --> mass mixing ratio --> volume mixing ratio
                q_ERA5_interp = (q_ERA5_interp / (1-q_ERA5_interp)) / epsilon # volume mixing ratio
   
                #Store in a Dataset
                
                ERA5_interp = xr.Dataset({'q': (['level'], q_ERA5_interp)},
                                         coords={'level': (['level'], level)})
                
                ERA5_interp["ta"]=(["level"], ta_ERA5_interp)


                # switch from altitude to pressure grid for sonde
                sonde = sonde.swap_dims({alt_var:p_var}).reset_coords()

                # Units conversion
                sonde[p_var] = sonde[p_var] * hPaToPa                       # convert from hPa to Pa 
                sonde[q_var] = (sonde[q_var] / (1- sonde[q_var])) / epsilon # calculate volume mixing ratio
                sonde[t_var] = sonde[t_var] + CtoK                          # C to K

                #in case indexes are not unique for the pressure (for radiosondes in particular), keep only the first index value
                _,index = np.unique(sonde[p_var], return_index=True)
                index = np.flip(index)
                sonde = sonde.isel({p_var:index})

                
                # Construct pressure grid: coarse where we have only the background sounding, changing abruptly to
                #   the a grid from min to max pressure of the sondes.

                # minimum pressure of sonde (in Pa)
                play_switch = np.ceil(sonde[p_var].min())

                # Choose background pressure from ERA grid. 
                # Then interpolate tropical atm profiles onto ERA grid
                # ERA in hPa, sonde/play_switch in Pa, convert ERA to Pa

                # background ERA pressure levels in Pa
                back_plays_ERA_Pa = q_ERA5['level'] * hPaToPa 
                back_plays_ERA  = np.sort(back_plays_ERA_Pa.where(back_plays_ERA_Pa < play_switch).dropna(dim='level'))
                # create 1hPa, 100Pa pressure grid for sonde
                sonde_plays = np.arange(play_switch, sonde[p_var].max(), deltaP)
                # append coarser ERA pressure grid to 1hPa sonde grid
                # play goes from low --> high pressure; top --> surface
                play = np.append(back_plays_ERA, sonde_plays) 

                # Interface pressures: mostly the average of the two neighboring layer pressures
                plev = np.append(np.append(back_plays_ERA.min(), 0.5 * (play[1:] + play[:-1])), play.max() + deltaP/2.)

                # Hydrostatic interpolation for the altitude
                
                top_level = sonde.where(sonde[p_var]==sonde[p_var].min(), drop=True)
                pres_top_sonde = top_level[p_var].values[0]
                alt_top_sonde = top_level[alt_var].values[0]
                temp_top_sonde = top_level[t_var].values[0]

                ERA5_interp = hydrostatic_interp(ERA5_interp, alt_top_sonde, pres_top_sonde, temp_top_sonde)
                ERA5_interp["SST"] = SST_ERA5_interp[0]


                zlay = np.append(ERA5_interp["zlay"].values, sonde[alt_var].interp({p_var:sonde_plays}))
                zlev = np.append(np.append(ERA5_interp["zlev"].max(), 0.5*(zlay[1:] + zlay[:-1])), zlay.min()/2)

                #Calculate solar angle using lat, lon, time
                alt_sol = pysolar.get_altitude(lat,lon, date_UTC)
                conv_deg_rad = np.pi/180
                cos_sza = np.sin(alt_sol*conv_deg_rad)  
                if(cos_sza <= 0): 
                    mu0 = 0.
                else:
                    mu0 = cos_sza

                    
                sfc_t = ERA5_interp.SST.values
                temp = np.append(ERA5_interp.ta.values, sonde[t_var].interp({p_var:sonde_plays}))            
                h2o = np.append(ERA5_interp.q.values, sonde[q_var].interp({p_var:sonde_plays}))

                profile = xr.Dataset({"launch_time":([], sonde.launch_time),\
                                      "platform":([], sonde.Platform.values),
                                      "tlay"   :(["play"], temp), \
                                      "play"   :(["play"], play), \
                                      "h2o":(["play"], h2o),  \
                                      "zlay":(["play"], zlay), \
                                      "plev"   :(["plev"], plev), \
                                      "zlev"   :(["plev"], zlev),\
                                      "sfc_emis":([], sfc_emis),\
                                      "sfc_alb":([], sfc_alb ),\
                                      "sfc_t":([], sfc_t),\
                                      "cos_sza":([], mu0),\
                                      "lw_dn"  :(["plev"], np.repeat(np.nan, plev.size)),\
                                      "lw_up"  :(["plev"], np.repeat(np.nan, plev.size)),\
                                      "lw_net" :(["plev"], np.repeat(np.nan, plev.size)),\
                                      "sw_dn"  :(["plev"], np.repeat(np.nan, plev.size)),\
                                      "sw_up"  :(["plev"], np.repeat(np.nan, plev.size)),\
                                      "sw_net" :(["plev"], np.repeat(np.nan, plev.size))})
#                 #
                # Add the other greenhouse gases
                # AL, haven't interpolated this onto ERA grid

                lowest = back_tropical_atm.p_lay.argmax()
                back_on_p = back_tropical_atm.swap_dims({'lay':'p_lay'}).reset_coords() # Background sounding on pressure layers
                back_on_p = back_on_p.rename({'p_lay':'play'}) # Rename p_lay into play
                for g in ghgs:
                    profile[g] = back_on_p["vmr_" + g].interp(play=play).fillna(back_tropical_atm["vmr_" + g].isel(lay=lowest))

                back_tropical_atm.close()
                sonde.close()

                all_profiles.append(profile)    
                
    return(all_profiles)

#functions to perform the interpolation

def interpScalar(values,lat_coord,lon_coord,lat,lon):    
    # We can think of changing option 'kind=' to something different than a linear interpolation
    return interp2d(lat_coord,lon_coord,values)(lat,lon)

def interpProfile(values,lat_coord,lon_coord,lat,lon):

    Nz = values.shape[0]
    profile_out = np.nan*np.zeros((Nz,))
    
    for i_z in range(Nz):
        profile_out[i_z] = interpScalar(values[i_z],lat_coord,lon_coord,lat,lon)
    
    return profile_out

def hydrostatic_interp(ERA5_interp, alt_top_sonde, pres_top_sonde, temp_top_sonde):
     
     #In ERA profile, we keep only pressure values inferior to the highest pressure measured by the sounding
    ERA5_interp = ERA5_interp.where(ERA5_interp.level<pres_top_sonde, drop=True)
    length = len(ERA5_interp.level)

    #array where we will store the altitude
    zlay = []
    zlev = []
        
    #We initialize this array
    first_level_above = ERA5_interp.isel(level=length-1)

    first_alt_above = alt_top_sonde + Rd/g*((first_level_above.ta.values + temp_top_sonde)/2\
                                            *math.log(pres_top_sonde/first_level_above.level.values))
    zlev.insert(0, first_alt_above)
    for i in range(length-1, 0, -1):
        profile_1 = ERA5_interp.isel(level=i)
        profile_2 = ERA5_interp.isel(level=i-1)

        new_alt = zlev[0] + Rd/g*((profile_1.ta.values + profile_2.ta.values)/2
                    *math.log(profile_1.level.values/profile_2.level.values))

        zlev.insert(0, new_alt)
        zlay.insert(0, (zlev[0]+zlev[1])/2)

    #for the top level, we assume same difference with the one before 
    zlay.insert(0, zlev[0])
    
    zlev = np.array(zlev)
    zlay = np.array(zlay)

    ERA5_interp["zlev"] = (["level"], zlev) 
    ERA5_interp["zlay"] = (["level"], zlay) 

    return ERA5_interp



# argparse variables: sonde_file, background_sounding_file, delta_p
# values when run in repo root
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Builds a netCDF file suitable for computing radiative fluxes by merging a background sounding a sonde file from Aspen")
    parser.add_argument("--sonde_file", type=str, default="input/HALO/20200119/AVAPS_Dropsondes/processed/D20200119_161410_PQC.nc",
                        help="Name of sonde file")
    parser.add_argument("--background_file", type=str, default='../input/tropical-atmosphere.nc',
                        help="Directory where reference values are")
    parser.add_argument("--ERA_dir", type=str, default='../input/',
                        help="Directory where ERA files are")
    parser.add_argument("--deltaP", type=int, default=100,
                        help="Pressure discretization of sonde (Pa, integer)")
    parser.add_argument("--sfc_emissivity", type=float, default=0.98, dest="emis",
                        help="Surface emissivity (spectrally constant)")
    parser.add_argument("--sfc_albedo", type=float, default=0.07, dest="alb",
                        help="Surface albedo (spectrally constant, same for direct and diffuse)")
    parser.add_argument("--cos_sza", type=float, default=0, dest="mu0",
                        help="Cosine of solar zenith angle, default is to compute from sonde file (someday)")
    parser.add_argument("--out_dir", type=str, default="../output/rad_profiles",
                        help="Directory where the output files should be saved")
    parser.add_argument("--out_file", type=str, default=None,
                        help="Output file name")
    args = parser.parse_args()

    # Generalize this
    output_dir  = args.out_dir
    
    all_profiles = combine_sonde_and_background(args.sonde_file, args.background_file, args.ERA_dir,\
                                           deltaP=args.deltaP, sfc_emis=args.emis, sfc_alb=args.alb, mu0=args.mu0)
    
    i=0
    for profile in all_profiles:
        #  str_launch_time = str(profile.launch_time.values)[:-10]
        Platform = str(profile.platform.values)
        output_file = "rrtmgp_" + '{:>04}'.format(i)+ "_%s.nc"%Platform
        print(output_file)
        profile.to_netcdf(os.path.join(output_dir, output_file))
        i+=1