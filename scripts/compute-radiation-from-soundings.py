#!/usr/env/python
import argparse
import xarray as xr
import numpy as np
import os

#
# Constants - could use a module if that's better
#
PaTohPa = 100.
epsilon = 0.6223 # Ratio of molar mass of water to dry air
CtoK    = 273.15 # Celsius to Kelvin
gtokg   = 1.e-3

def combine_sonde_and_background(sonde_file, background_file, deltaP=100, sfc_emis=.98, sfc_alb=0.07,mu0=1.):
    #
    # Maybe someone can show me how to use Python doc strings?
    #   This routine takes a single background sounding extracted from the Garand atmospheres distributed with RTE+RRTMGP
    #   and pastes it above a sounding file following the "QC netcdf" conventions from the Apsen software
    #
    #
    # Background sounding
    #
    back = xr.open_dataset(background_file)
    #
    # Parse the sondes, dropping levels at which any of
    #    temperature, water vapor mass mixing ratio, pressure, or altitude are missing
    #
    sonde = xr.open_dataset(sonde_file).dropna(dim='time',subset=['time']). \
                                        swap_dims({'time':'pres'}).reset_coords(). \
                                        dropna(dim='pres',subset=['alt','pres','mr','tdry'],how='any')
    #
    # Units conversion
    #
    sonde["pres"] = sonde.pres * PaTohPa       # hPa to Pa
    sonde["mr"]   = sonde.mr * gtokg / epsilon # Mass to volume mixing ratio
    sonde["tdry"] = sonde.tdry + CtoK          # C to K

    #
    # Construct pressure grid: coarse where we have only the background sounding, changing abruptly to
    #   the a grid from min to max pressure of the sondes.
    #
    play_switch = np.ceil(sonde.pres.min())
    # Layer pressures from background sounding in increasing order
    back_plays  = np.sort(back.p_lay.where(back.p_lay < play_switch).dropna(dim='lay'))
    sonde_plays = np.arange(play_switch, sonde.pres.max(), deltaP)
    nback       = back_plays.size
    play = np.append(back_plays, sonde_plays)
    #
    # Interface pressures: mostly the average of the two neighboring layer pressures
    #
    plev = np.append(np.append(back.p_lev.min(), 0.5 * (play[1:] + play[:-1])), play.max() + deltaP/2.)

    #
    # Interpolate values onto new grid
    #
    temp = np.append(back.swap_dims({'lay':'p_lay'}).reset_coords().t_lay.interp(p_lay=back_plays), \
                     sonde.tdry.interp(pres=sonde_plays))
    h2o  = np.append(back.swap_dims({'lay':'p_lay'}).reset_coords().vmr_h2o.interp(p_lay=back_plays), \
                     sonde.mr.interp(pres=sonde_plays))
    #
    # Index with the greatest pressure - where pressures in the sonde are higher than any in the background, use the
    #   value from the highest pressure/lowest level
    #
    lowest = back.p_lay.argmax()
    o3   = back.swap_dims({'lay':'p_lay'}).reset_coords().vmr_o3.interp(p_lay=play).fillna(back.vmr_o3.isel(lay=lowest))
    n2   = back.swap_dims({'lay':'p_lay'}).reset_coords().vmr_n2.interp(p_lay=play).fillna(back.vmr_n2.isel(lay=lowest))
    o2   = back.swap_dims({'lay':'p_lay'}).reset_coords().vmr_o2.interp(p_lay=play).fillna(back.vmr_o2.isel(lay=lowest))
    co   = back.swap_dims({'lay':'p_lay'}).reset_coords().vmr_co.interp(p_lay=play).fillna(back.vmr_co.isel(lay=lowest))
    co2  = back.swap_dims({'lay':'p_lay'}).reset_coords().vmr_co2.interp(p_lay=play).fillna(back.vmr_co2.isel(lay=lowest))
    ch4  = back.swap_dims({'lay':'p_lay'}).reset_coords().vmr_ch4.interp(p_lay=play).fillna(back.vmr_ch4.isel(lay=lowest))
    n2o  = back.swap_dims({'lay':'p_lay'}).reset_coords().vmr_n2o.interp(p_lay=play).fillna(back.vmr_n2o.isel(lay=lowest))

    if(mu0 <= 0): mu0 = 1. # replace with computation from sonde date, time, lat/lon

    profile = xr.Dataset({"tlay"   :(["play"], temp), \
                          "play"   :(["play"], play), \
                          "vmr_h2o":(["play"], h2o),  \
                          "vmr_co2":(["play"], co2),  \
                          "vmr_ch4":(["play"], ch4),  \
                          "vmr_n2o":(["play"], n2o),  \
                          "vmr_o2" :(["play"], o2 ),  \
                          "vmr_o3" :(["play"], o3 ),  \
                          "vmr_n2" :(["play"], n2 ),  \
                          "vmr_co" :(["play"], co),   \
                          "plev"   :(["plev"], plev), \
                          "sfc_emis":([], sfc_emis),  \
                          "sfc_alb":([], sfc_alb ),  \
                          "cos_sza":([], mu0),       \
                          "lw_dn"  :(["plev"], np.repeat(np.nan, plev.size)),\
                          "lw_up"  :(["plev"], np.repeat(np.nan, plev.size)),\
                          "lw_net" :(["plev"], np.repeat(np.nan, plev.size)),\
                          "sw_dn"  :(["plev"], np.repeat(np.nan, plev.size)),\
                          "sw_up"  :(["plev"], np.repeat(np.nan, plev.size)),\
                          "sw_net" :(["plev"], np.repeat(np.nan, plev.size))})
    back.close()
    sonde.close()
    return(profile)

# argparse variables: sonde_file, background_sounding_file, delta_p
# values when run in repo root
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Builds a netCDF file suitable for computing radiative fluxes by merging a background sounding a sonde file from Aspen")
    parser.add_argument("--sonde_file", type=str, default="input/HALO/20200119/AVAPS_Dropsondes/processed/D20200119_161410_PQC.nc",
                        help="Name of sonde file")
    parser.add_argument("--background_file", type=str, default='tropical-atmosphere.nc',
                        help="Directory where reference values are")
    parser.add_argument("--deltaP", type=int, default=100,
                        help="Pressure discretization of sonde (Pa, integer)")
    parser.add_argument("--sfc_emissivity", type=float, default=0.98, dest="emis",
                        help="Surface emissivity (spectrally constant)")
    parser.add_argument("--sfc_albedo", type=float, default=0.07, dest="alb",
                        help="Surface albedo (spectrally constant, same for direct and diffuse)")
    parser.add_argument("--cos_sza", type=float, default=0, dest="mu0",
                        help="Cosine of solar zenith angle, default is to compute from sonde file (someday)")
    args = parser.parse_args()

    # Generalize this
    output_dir      = '.'
    output_file     = os.path.basename(args.sonde_file)[:-3] + "_rad.nc"

    # Any error checking on arguments? 

    profile = combine_sonde_and_background(args.sonde_file, args.background_file, \
                                           deltaP=args.deltaP, sfc_emis=args.emis, sfc_alb=args.alb, mu0=args.mu0)
    profile.to_netcdf(os.path.join(output_dir, output_file))
