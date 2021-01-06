#!/bin/bash

# by Ben Fildier and Ludovic Touze-Peiffer
# Uses radiative transfer code by Robert Pincus to calculate
# radiative profiles from HALO dropsonde data


wdir=${PWD%/*}
cdir=${wdir}/code

# Replace with the path to your input file   

#---- Ludo

#ifile=/media/ludo/DATA/google-drive/Thèse/EUREC4a/github/Input/all_sondes_2.nc
#odir=${wdir}/output/rad_profiles-test
#ERAdir=/media/ludo/DATA/google-drive/Thèse/EUREC4a/github/Input/ERA
#ifile=/run/media/ludo/DATA/google-drive/Thèse/EUREC4a/github/Input/Dropsondes/all_dropsondes.nc
#odir=${wdir}/output/rad_profiles_dropsondes

#--- Ben

# directory structure for perturbations
idir=/Users/bfildier/Code/analyses/EUREC4A/EUREC4A_rad_profiles/input/perturbations
odir=/Users/bfildier/Code/analyses/EUREC4A/EUREC4A_rad_profiles/output/perturbations
#subdirs="baseline ERA_SST_m021K ERA_SST_p021K ERA_q_m30pct ERA_q_p30pct sondes_q_m3pct sondes_q_p3pct"
subdirs="ERA_SST_m042K ERA_SST_p042K"

#ERAdir=/Users/bfildier/Data/EUREC4A/ERA
#ifile=/Users/bfildier/Data/EUREC4A/merged/sondes/all_sondes.nc
#ifile=/Users/bfildier/Data/EUREC4A/merged/sondes/all_sondes_m3pct.nc
#odir=${wdir}/output/rad_profiles-test

for subdir in `echo $subdirs`; do

	mkdir -p ${odir}

	# python combine_ERA_and_sonde_profiles.py --sonde_file=${ifile} --out_dir=${odir} --ERA_dir=${ERAdir}
#	python combine_ERA_and_sonde_profiles.py --sonde_file=${idir}/${subdir}/all_sondes.nc --out_dir=${odir}/${subdir} --ERA_dir=${idir}/${subdir}

#	for ofile in `ls ${odir}/${subdir}/rrtmgp_0???_*.nc`; do	
#		echo 'Compute radiation profile '$ofile
#             	${cdir}/sonde_radiation $ofile
#             	echo " "
#	done

#	for ofile in `ls ${odir}/${subdir}/rrtmgp_1???_*.nc`; do
#	        echo 'Compute radiation profile '$ofile
#             	${cdir}/sonde_radiation $ofile
#            	echo " "
#	done

#	for ofile in `ls ${odir}/${subdir}/rrtmgp_2???_*.nc`; do
#             	echo 'Compute radiation profile '$ofile
#             	${cdir}/sonde_radiation $ofile
#             	echo " "
#	done

	#python post_processing.py --in_dir=${odir} --out_dir=$odir --comp_qrad=True
	python post_processing.py --in_dir=${odir}/${subdir} --out_dir=${odir}/${subdir} --comp_qrad=True

#To compute the radiative heating itself, run instead

#python post_processing.py --in_dir="../output/rad_profiles" --out_dir="../output/" --comp_qrad=True

done

exit 0
