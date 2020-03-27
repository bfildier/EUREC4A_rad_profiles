#!/bin/bash

# by Ben Fildier and Ludovic Touze-Peiffer
# Uses radiative transfer code by Robert Pincus to calculate
# radiative profiles from HALO dropsonde data


wdir=${PWD%/*}
cdir=${wdir}/code

# Replace with the path to your input file   

#ifile=/run/media/ludo/DATA/google-drive/Thèse/EUREC4a/github/Input/Radiosondes/all_radiosondes.nc
#odir=${wdir}/output/rad_profiles_radiosondes

ifile=/run/media/ludo/DATA/google-drive/Thèse/EUREC4a/github/Input/Dropsondes/all_dropsondes.nc
odir=${wdir}/output/rad_profiles_dropsondes

mkdir -p ${odir}

python combine_reference_and_sonde_profiles.py --sonde_file=${ifile} --out_dir=${odir}

for ofile in `ls ${odir}/*.nc`; do
               echo 'Compute radiation profile '$ofile
               ${cdir}/sonde_radiation $ofile
               echo " "
done


python post_processing.py --in_dir=${odir} --out_dir="../output/" --comp_qrad=False

#To compute the radiative heating itself, run instead

#python post_processing.py --in_dir="../output/rad_profiles" --out_dir="../output/" --comp_qrad=True

