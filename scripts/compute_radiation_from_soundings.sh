#!/bin/bash

# by Ben Fildier and Ludovic Touze-Peiffer
# Uses radiative transfer code by Robert Pincus to calculate
# radiative profiles from HALO dropsonde data


wdir=${PWD%/*}
cdir=${wdir}/code

# Replace with the path to your input file   

#---- Ludo

ifile=/run/media/ludo/DATA/google-drive/Thèse/EUREC4a/github/Input/all_sondes_2.nc
odir=${wdir}/output/rad_profiles
#ifile=/run/media/ludo/DATA/google-drive/Thèse/EUREC4a/github/Input/Dropsondes/all_dropsondes.nc
#odir=${wdir}/output/rad_profiles_dropsondes

#--- Ben
#ifile=/Users/bfildier/Data/EUREC4A/merged/sondes/all_radiosondes.nc
#odir=${wdir}/output/rad_profiles_radiosondes
#ifile=/Users/bfildier/Data/EUREC4A/merged/sondes/all_dropsondes.nc
#odir=${wdir}/output/rad_profiles_dropsondes
#ifile=/Users/bfildier/Data/EUREC4A/merged/sondes/proxy_sondes_fixedT.nc
#odir=${wdir}/output/rad_profiles

mkdir -p ${odir}

python combine_reference_and_sonde_profiles.py --sonde_file=${ifile} --out_dir=${odir}

#for ofile in `ls ${odir}/*.nc`; do
#              echo 'Compute radiation profile '$ofile
#              ${cdir}/sonde_radiation $ofile
#              echo " "
#done
#
#python post_processing.py --in_dir=${odir} --out_dir=$odir --comp_qrad=True

#To compute the radiative heating itself, run instead

#python post_processing.py --in_dir="../output/rad_profiles" --out_dir="../output/" --comp_qrad=True

