#!/usr/env/python
from pysolar.solar import *
import datetime
import argparse
import pytz
import xarray as xr
import numpy as np
import os
import pytz
from metpy import calc as mpcalc


#
# Constants - could use a module if that's better
#
hPatoPa = 100.
epsilon = 0.6223 # Ratio of molar mass of water to dry air
CtoK    = 273.15 # Celsius to Kelvin
gtokg   = 1.e-3
ghgs    = ["co2", "ch4", "n2o", "o3", "o2", "n2", "co"]
deltaZ = 8.62                     #in m corresponding to deltaP=100Pa for T=25°C 
#
# Need to add surface temperature for LW
#

def combine_sonde_and_background(all_sondes_file, background_file, SST_dir, deltaP=100, sfc_emis=.98, sfc_alb=0.07, mu0=1., ghgs=ghgs):
    #
    # Maybe someone can show me how to use Python doc strings?
    #   This routine takes a single background sounding extracted from the Garand atmospheres distributed with RTE+RRTMGP
    #   and pastes it above a sounding file following the "QC netcdf" conventions from the Apsen software
    #
    #
    # Background sounding
    #
    back = xr.open_dataset(background_file)
    
    back["p_lay"].attrs['units'] = 'Pa'
    
    back_zlay = mpcalc.pressure_to_height_std(back["p_lay"])
    back["zlay"] = (["lay"] ,back_zlay.magnitude*1000)
        
    back_zlev = mpcalc.pressure_to_height_std(back["p_lev"])
    back["zlev"] = (["lev"] ,back_zlev.magnitude*1000)
    #
    # Parse the sonde, dropping levels at which any of
    #    temperature, water vapor mass mixing ratio, pressure, or altitude are missing
    #    and making pressure the coordinate
    #
    
    all_sondes = xr.open_dataset(all_sondes_file)
    all_profiles = []
    
    number_sondes = len(all_sondes.launch_time)
    
    for i in range(number_sondes):
   # for i in range(100,106):

        alt_var = "height"
        p_var = 'pressure'
        q_var = 'specific_humidity'
        t_var = 'temperature'
        rh_var = 'relative_humidity'

        sonde = all_sondes.isel(launch_time = i).dropna(dim=alt_var,\
                               subset=[alt_var, p_var,q_var,t_var],\
                                                   how="any")
        
        if (sonde[alt_var].values.size < 10 or sonde[alt_var].min() > 200):
            print("The sonde is empty ")
            sonde.close()            
            
        else:
            
            
            #We get rid of z=0 due to interface issue
          #  sonde = sonde.where(sonde.gpsalt > 0, drop = True)
        
            sonde = sonde.swap_dims({alt_var:p_var}).reset_coords()
        #
        # Units conversion
        #
            sonde[p_var] = sonde[p_var] * hPatoPa       # hPa to Pa
            # sonde[q_var] = sonde[q_var] * gtokg / epsilon # Mass to volume mixing ratio
            sonde[q_var] = sonde[q_var]/(1-sonde[q_var]) / epsilon # specific humidity to volume mixing ratio
            sonde[t_var] = sonde[t_var] + CtoK          # C to K
    
    
            #only for radiosondes, otherwise comment following lines
    
            _,index = np.unique(sonde[p_var], return_index=True)

            index = np.flip(index)
       
            sonde = sonde.isel({p_var:index})
                  
            #
            # Construct pressure grid: coarse where we have only the background sounding, changing abruptly to
            #   the a grid from min to max pressure of the sondes.
            #
            play_switch = np.ceil(sonde[p_var].min())
            # Layer pressures from background sounding in increasing order
            back_plays  = np.sort(back.p_lay.where(back.p_lay < play_switch).dropna(dim='lay'))
            sonde_plays = np.arange(play_switch, sonde[p_var].max(), deltaP)
            play = np.append(back_plays, sonde_plays)
            
            #
            # Interface pressures: mostly the average of the two neighboring layer pressures
            #
            plev = np.append(np.append(back.p_lev.min(), 0.5 * (play[1:] + play[:-1])), play.max() + deltaP/2.)

            #
            # Interpolate values onto new grid
            #
            
            back_zlay = back.swap_dims({'lay':'p_lay'}).reset_coords().zlay.interp(p_lay=back_plays)
            back_zlay = back_zlay + sonde[alt_var].max().values - back_zlay.min().values + 100
            
            zlay = np.append(back_zlay,sonde[alt_var].interp({p_var:sonde_plays}))
            zlev = np.append(np.append(back.zlev.max(), 0.5*(zlay[1:] + zlay[:-1])), zlay.min() - deltaZ/2)
                                    
            temp = np.append(back.swap_dims({'lay':'p_lay'}).reset_coords().t_lay.interp(p_lay=back_plays), \
                         sonde[t_var].interp({p_var:sonde_plays}))
            
            
            h2o  = np.append(back.swap_dims({'lay':'p_lay'}).reset_coords().vmr_h2o.interp(p_lay=back_plays), \
                         sonde[q_var].interp({p_var:sonde_plays}))

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
            
                date = datetime.datetime.strptime(str(sonde.launch_time.values).split('.')[0],'%Y-%m-%dT%H:%M:%S')
                timezone = pytz.timezone("UTC")
                date = timezone.localize(date)
                alt_sol = get_altitude(lat,lon, date)
                conv_deg_rad = np.pi/180
                cos_sza = np.sin(alt_sol*conv_deg_rad)  

                if(cos_sza <= 0): 
                    mu0 = 0.
                else:
                    mu0 = cos_sza

                date_strg_SST = date.strftime("%Y-%m-%d")
                SST_path = os.path.join(SST_dir, "dataset-catsat-nrt-global-infrared-sst-hr-daily-" + date_strg_SST + ".nc")

                SST_file = xr.open_dataset(SST_path)
                SST_file = SST_file["Grid_0001"]
                SST_file = SST_file.interpolate_na(dim="NbLongitudes")
                SST_file = SST_file.interpolate_na(dim="NbLatitudes")

                SST_file = SST_file.sel(NbLatitudes = lat.values, method="nearest", drop=True).dropna(dim="NbLongitudes")
                SST_file = SST_file.sel(NbLongitudes = lon.values, method="nearest",drop=True)

                sfc_t = SST_file.values[0] + CtoK

           # replace with computation from sonde date, time, lat/lon

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
                #
                lowest = back.p_lay.argmax()
                back_on_p = back.swap_dims({'lay':'p_lay'}).reset_coords() # Background sounding on pressure layers
                back_on_p = back_on_p.rename({'p_lay':'play'}) # Rename p_lay into play
                for g in ghgs:
                    profile[g] = back_on_p["vmr_" + g].interp(play=play).fillna(back["vmr_" + g].isel(lay=lowest))

                back.close()
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
