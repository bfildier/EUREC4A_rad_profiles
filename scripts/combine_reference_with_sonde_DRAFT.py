#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

# =============================================================================
# To do's
# =============================================================================
# add minimum alt for sonde to arg parse? Line 92

# use hydrostatic balance to calculate height from top sonde alt. value --> ERA5 alt top 

#SST_sel_time = SST.sel(time=date_local, method="nearest")
                # interpolate around local lat/lon point
#%%
#!/usr/env/python
import pysolar.solar as pysolar
import datetime
import argparse
import pytz
import xarray as xr
import numpy as np
import os
import pytz
from metpy import calc as mpcalc

#%%
# load ERA5 SST, air temperature, and specific humidity hourly fields

# AL, started adding these files to Arg parse, but just load locally for now

input_dir = '/Users/annaleaalbright/Dropbox/EUREC4A/RadiativeProfiles/Data/'
q_ERA5 = xr.open_dataset(input_dir + 'q_2020_01_02_ERA5_hourly_Barbados.nc')
ta_ERA5 = xr.open_dataset(input_dir + 'ta_2020_01_02_ERA5_hourly_Barbados.nc')
SST_ERA5 = xr.open_dataset(input_dir + 'SST_2020_01_02_ERA5_hourly_Barbados.nc')

#%%

# Define parameters

hPaToPa = 100.                                          # Convert from hPa to Pa
PatohPa = 1 / 100                                       # Convert from Pa to hPa
epsilon = 0.6223                                        # Ratio of molar mass of water to dry air
CtoK    = 273.15                                        # Celsius to Kelvin
gtokg   = 1.e-3                                         # g/kg to kg/kg
ghgs    = ["co2", "ch4", "n2o", "o3", "o2", "n2", "co"] # greenhouse gases
deltaZ = 8.62                                           #in m corresponding to deltaP=100Pa for T=25°C 


def combine_sonde_and_background(all_sondes_file, background_file, deltaP=100, sfc_emis=.98, sfc_alb=0.07, mu0=1., ghgs=ghgs):

    # Background sounding
    back_tropical_atm = xr.open_dataset(background_file)
    
    # add filenames to Argparse
    ta_ERA5 = xr.open_dataset(ta_ERA5_file)
    SST_ERA5 = xr.open_dataset(SST_ERA5_file)
    q_ERA5 = xr.open_dataset(q_ERA5_file)
    # convert from specific humidity to volume mixing ratio
    # AL, does someone want to double check this equation, just given its importance? :D
    #convert specific humidity --> mass mixing ratio --> volume mixing ratio
    vmr_ERA5 = (q_ERA5 / (1-q_ERA5)) / epsilon # volume mixing ratio

    
    
    # Ludo, use hydrostatic balance to calculate height from top sonde alt. value --> ERA5 alt top 


    
    #
    # Parse the sonde, dropping levels at which any of
    #    temperature, water vapor mass mixing ratio, pressure, or altitude are missing
    #    and making pressure the coordinate
    
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
        min_alt = 3000 # add to arg parse?
        if (sonde[alt_var].values.size < 10 or sonde[alt_var].min() > 200 or sonde[alt_var].max() < min_alt ):
            print("The sonde is empty ")
            sonde.close()            
            
        else:
            
            
            #We get rid of z=0 due to interface issue
          #  sonde = sonde.where(sonde.gpsalt > 0, drop = True)
             
           # switch from altitude to pressure grid for sonde
            sonde = sonde.swap_dims({alt_var:p_var}).reset_coords()

            # Units conversion
            # AL, double check, conversion factor was named PatohPa but number was right to go from hPa--> Pa
            # renamed hPatoPa
            
            sonde[p_var] = sonde[p_var] * hPaToPa                       # convert from hPa to Pa 
            #sonde[p_var].attrs['units'] = 'Pa'
            sonde[q_var] = (sonde[q_var] / (1- sonde[q_var])) / epsilon # calculate volume mixing ratio
                                                                        # specific humidity q --> mass mixing ratio m: m = q/(1-q)
                                                                        # VMR = MMR / epsilon, where epsilon = Rw/Rd
            #sonde[q_var].attrs['long_name'] = 'volume mixing ratio'
            sonde[t_var] = sonde[t_var] + CtoK                          # C to K
            #sonde[t_var].attrs['units'] = 'K'
    
            #only for radiosondes, otherwise comment following lines
    
#            _,index = np.unique(sonde[p_var], return_index=True)
#
#            index = np.flip(index)
#       
#            sonde = sonde.isel({p_var:index})
                  
            #
            # Construct pressure grid: coarse where we have only the background sounding, changing abruptly to
            #   the a grid from min to max pressure of the sondes.
            
            # minimum pressure of sonde (in Pa)
            play_switch = np.ceil(sonde[p_var].min())
            
            # Choose background pressure from ERA grid. 
            # Then interpolate tropical atm profiles onto ERA grid
            # ERA in hPa, sonde/play_switch in Pa, convert ERA to Pa
    
            #back_plays  = np.sort(back.p_lay.where(back.p_lay < play_switch).dropna(dim='lay'))
            back_plays_ERA_Pa = q_ERA5['level'] * hPaToPa 
            # background ERA pressure levels in Pa
            back_plays_ERA  = np.sort(back_plays_ERA_Pa.where(back_plays_ERA_Pa < play_switch).dropna(dim='level'))
            # create 1hPa, 100Pa pressure grid for sonde
            sonde_plays = np.arange(play_switch, sonde[p_var].max(), deltaP)
            # append coarser ERA pressure grid to 1hPa sonde grid
            # play goes from low --> high pressure; top --> surface
            play = np.append(back_plays_ERA, sonde_plays) 
            
            # Interface pressures: mostly the average of the two neighboring layer pressures
            plev = np.append(np.append(back_plays_ERA.min(), 0.5 * (play[1:] + play[:-1])), play.max() + deltaP/2.)


            # Interpolate values onto new altitude grid

            
            #back_zlay = back.swap_dims({'lay':'p_lay'}).reset_coords().zlay.interp(p_lay=back_plays)
            #back_zlay = back_zlay + sonde[alt_var].max().values - back_zlay.min().values + 100
            
            # AL, got an error message here
            # "InvalidIndexError: Reindexing only valid with uniquely valued Index objects'
            # Ludo will add height from hydrostatic balance
            zlay = np.append(back_zlay,sonde[alt_var].interp({p_var:sonde_plays}))
            zlev = np.append(np.append(back.zlev.max(), 0.5*(zlay[1:] + zlay[:-1])), zlay.min() - deltaZ/2)

        #
        # Index with the greatest pressure - where pressures in the sonde are higher than any in the background, use the
        #   value from the highest pressure/lowest level
        #
            lat = sonde.latitude.dropna(dim=p_var)[0]
            lon = sonde.longitude.dropna(dim=p_var)[0]
            
            if (lat > 20 or lat < 4 or lon < -65 or lon > -50):
                print("outside SST bounds")
                sonde.close()
                
            else:
                fmt = '%Y-%m-%dT%H:%M:%S'
                date = datetime.datetime.strptime(str(sonde.launch_time.values).split('.')[0], fmt)
                date_UTC = pytz.timezone('UTC').localize(date)
                #date_local = date_UTC.astimezone(pytz.timezone('America/Barbados'))  # convert from UTC to local time

                alt_sol = pysolar.get_altitude(lat,lon, date_UTC)
                conv_deg_rad = np.pi/180
                cos_sza = np.sin(alt_sol*conv_deg_rad)  
                if(cos_sza <= 0): 
                    mu0 = 0.
                else:
                    mu0 = cos_sza
                        
                # =============================================================================
                #                  select SST value
                # =============================================================================

                # ERA5 data in UTC
                SST_sel_time = SST_ERA5.sel(time=date_UTC, method="nearest")
                
                # interpolate around local lat/lon point
                
                # find local SST value
                #SST_sel_lat = SST_sel_time.sel(latitude = lat.values, method="nearest", drop=True).dropna(dim="longitude")
                #SST_sel_lat_lon = SST_sel_lat.sel(longitude = lon.values, method="nearest",drop=True)
                
                sfc_t = SST_sel_lat_lon.sstk.values # in K

                # =============================================================================
                #           select ERA5 air temperature profile 
                #           and append to sonde
                # =============================================================================
                # (level 37, latitude 45, longitude 45) 
                ta_sel_time = ta_ERA5.sel(time=date_UTC, method="nearest")

                # 2d interpolation around local lat/lon point
                
                #ta_sel_lat = ta_sel_time.sel(latitude = lat.values, method="nearest", drop=True).dropna(dim="longitude")
                #ta_sel_lat_lon = ta_sel_lat.sel(longitude = lon.values, method="nearest",drop=True)
                # level = pressure level, air_pressure hPa -- from 1-->1000 hPa
                # array: ta_sel_lat_lon.ta.values
                
                # append air temperature from ERA5 to sonde temperature 
                # AL, back_plays_ERA and sonde_plays both in Pa
                # choose interpolated temperature, don't need to switch dims because 
                # it should already be a pressure grid
                temp = np.append(ta_ERA_interp, sonde[t_var].interp({p_var:sonde_plays}))

                #temp = np.append(back.swap_dims({'lay':'p_lay'}).reset_coords().t_lay.interp(p_lay=back_plays), \
                #         sonde[t_var].interp({p_var:sonde_plays}))
            

                # =============================================================================
                #             select ERA5 volume mixing ratio
                # =============================================================================
                
                vmr_sel_time = vmr_ERA5.sel(time=date_UTC, method="nearest")
                
                # 2d interpolation, select lat/lon point, append to sonde
                
                h2o = np.append(vmr_ERA_interp, sonde[q_var].interp({p_var:sonde_plays}))

                #h2o  = np.append(back.swap_dims({'lay':'p_lay'}).reset_coords().vmr_h2o.interp(p_lay=back_plays), \
                         #sonde[q_var].interp({p_var:sonde_plays}))


                profile = xr.Dataset({"launch_time":([], sonde.launch_time),\
                                      "platform":([], sonde.Platform.values),
                                      "tlay"   :(["play"], temp), \
                                      "play"   :(["play"], play), \
                                      "h2o":(["play"], h2o),  \
                                      "zlay":(["play"], zlay),  \
                                      "plev"   :(["plev"], plev), \
                                      "zlev"   :(["plev"], zlev), \
                                      "sfc_emis":([], sfc_emis),  \
                                      "sfc_alb":([], sfc_alb ),  \
                                      "sfc_t":([], sfc_t),  \
                                      "cos_sza":([], mu0),  \
                                      "lw_dn"  :(["plev"], np.repeat(np.nan, plev.size)),\
                                      "lw_up"  :(["plev"], np.repeat(np.nan, plev.size)),\
                                      "lw_net" :(["plev"], np.repeat(np.nan, plev.size)),\
                                      "sw_dn"  :(["plev"], np.repeat(np.nan, plev.size)),\
                                      "sw_up"  :(["plev"], np.repeat(np.nan, plev.size)),\
                                      "sw_net" :(["plev"], np.repeat(np.nan, plev.size))})
                #
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

# argparse variables: sonde_file, background_sounding_file, delta_p
# values when run in repo root
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Builds a netCDF file suitable for computing radiative fluxes by merging a background sounding a sonde file from Aspen")
    parser.add_argument("--sonde_file", type=str, default="input/HALO/20200119/AVAPS_Dropsondes/processed/D20200119_161410_PQC.nc",
                        help="Name of sonde file")
    parser.add_argument("--background_file", type=str, default='../input/tropical-atmosphere.nc',
                        help="Directory where reference values are")
    
    # Add ERA5 files to input folder
    parser.add_argument("--q_ERA5", type=str, default='../input/q_2020_01_02_ERA5_hourly_Barbados.nc',
                        help="Name of ERA5 specific humidity file")
    parser.add_argument("--ta_ERA5", type=str, default='../input/ta_2020_01_02_ERA5_hourly_Barbados.nc',
                        help="Name of ERA5 air temperature file")
    parser.add_argument("--SST_ERA5", type=str, default='../input/SST_2020_01_02_ERA5_hourly_Barbados.nc',
                        help="Name of ERA5 SST file")

    parser.add_argument("--SST_dir", type=str, default='../input/SST/',
                        help="Directory where SST netcdf data are")
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
    
    all_profiles = combine_sonde_and_background(args.sonde_file, args.background_file, args.SST_dir,\
                                           deltaP=args.deltaP, sfc_emis=args.emis, sfc_alb=args.alb, mu0=args.mu0)
    
    i=0
    for profile in all_profiles:
        #  str_launch_time = str(profile.launch_time.values)[:-10]
        Platform = str(profile.platform.values)
        output_file = "rrtmgp_" + '{:>04}'.format(i)+ "_%s.nc"%Platform
        print(output_file)
        profile.to_netcdf(os.path.join(output_dir, output_file))
        i+=1