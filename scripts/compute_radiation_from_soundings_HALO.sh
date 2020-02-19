#!/bin/bash

# by Ben Fildier and Ludovic Touze-Peiffer
# Uses radiative transfer code by Robert Pincus to calculate
# radiative profiles from HALO dropsonde data


daylist=(20200119 20200124 20200126 20200128 20200130 20200131 20200202 20200205 20200207)

download='false'

# Replace for your own input directory 
idir=/Users/bfildier/Data/EUREC4A/aeris/eurec4a-data/AIRCRAFT
if [ "$download" == 'true' ]; then
	idir=${wdir}/input/Data_local/Measurements
fi

wdir=${PWD%/*}
odir=${wdir}/output
cdir=${wdir}/code
sdir=${wdir}/scripts

if [ "$download" == 'true' ]; then

	wget -r -nH --no-parent --cut-dirs=3 -P ${idir} https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/HALO/AVAP-DROPSONDES/2020/ 

fi

for iday in {0..8} ; do

	day=${daylist[iday]} 
	echo ${day}
	idir_day=${idir}/HALO/AVAP-DROPSONDES/${day:0:4}/RF??_${day}
	odir_day=${odir}/${day}

	cd $sdir

	# Creates file and combines sonde profile with standard atmosphere
	for ifile in `ls ${idir_day}/*PQC.nc ${idir_day}/processed/*PQC.nc`; do 

		filename=${ifile##*\/}	
		echo ${filename}
		echo 'Combining sonde with background '$ifile
		mkdir -p ${odir_day}

		python combine_reference_and_sonde_profiles.py --sonde_file=${ifile} --out_dir=${odir_day} 

		echo " "	
	done

	# Calculates radiative profiles
	for ofile in `ls ${odir_day}/*.nc`; do
		echo 'Compute radiation profile '$ofile
		${cdir}/sonde_radiation $ofile
		echo " "
	done
		
done



#ulimit -s unlimited
#sonde_path=${1%\/*}
#sonde_filename=${1##*\/}
#echo 'processing sonde '$sonde_filename
#rad_filename=${sonde_filename%*.nc}'_rad.nc'
#
## Merge dropsonde profiles with standard atmosphere profiles
#echo 'merging sounding and reference profiles in new file '$rad_filename
#
## Compute radiative profles
#echo "compute radiative profiles and append them to ${rad_filename}" 
#../code/sonde_radiation ${rad_filename}
#
##echo "move ${rad_filename} to ${sonde_path}"
#mv ${rad_filename} ${sonde_path}/
