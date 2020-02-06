#!/bin/bash

# by Ben Fildier and Ludovic Touze-Peiffer
# Uses radiative transfer code by Robert Pincus to calculate
# radiative profiles from HALO dropsonde data

sonde_path=${1%\/*}
sonde_filename=${1##*\/}
echo 'processing sonde '$sonde_filename
rad_filename=${sonde_filename%*.nc}'_rad.nc'

# Merge dropsonde profiles with standard atmosphere profiles
echo 'merging sounding and reference profiles in new file '$rad_filename
python combine_reference_and_sonde_profiles.py --sonde_file ${sonde_path}/${sonde_filename}

# Compute radiative profles
echo "compute radiative profiles and append them to ${rad_filename}" 
../code/sonde_radiation ${rad_filename}

#echo "move ${rad_filename} to ${sonde_path}"
#mv ${rad_filename} ${sonde_path}/
