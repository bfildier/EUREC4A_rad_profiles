#!/bin/bash

# True location of input data
ERAdir=/Users/bfildier/Data/EUREC4A/ERA/
sondedir=/Users/bfildier/Data/EUREC4A/merged/sondes/
currentdir=/Users/bfildier/Code/analyses/EUREC4A/EUREC4A_rad_profiles/scripts/
#subdirs='baseline ERA_SST_m021K ERA_SST_p021K ERA_q_m30pct ERA_q_p30pct sondes_q_m3pct sondes_q_p3pct'
subdirs='ERA_SST_m042K ERA_SST_p042K'


## Create directory tree
# in input folder
mkdir ${currentdir}/../input/perturbations
cd ${currentdir}/../input/perturbations
mkdir $subdirs
# in output folder
mkdir ${currentdir}/../output/perturbations
cd ${currentdir}/../output/perturbations
mkdir $subdirs

## Create symbolic links to input data in each case
# link default files first
for subdir in `echo $subdirs`; do

	cd ${currentdir}/../input/perturbations/$subdir
	# sondes
	ln -s ${sondedir}/all_sondes.nc all_sondes.nc
	# ERA data
	for varid in `echo SST q ta`; do
 		ln -s ${ERAdir}/${varid}_2020_01_02_ERA5_hourly_Barbados.nc ${varid}_2020_01_02_ERA5_hourly_Barbados.nc
	done
done

## link perturbed files for relevant directories

# SST lower
cd ${currentdir}/../input/perturbations/ERA_SST_m021K
rm SST_2020_01_02_ERA5_hourly_Barbados.nc
ln -s ${ERAdir}/SST_2020_01_02_ERA5_hourly_Barbados_m021K.nc SST_2020_01_02_ERA5_hourly_Barbados.nc

# SST upper
cd ${currentdir}/../input/perturbations/ERA_SST_p021K
rm SST_2020_01_02_ERA5_hourly_Barbados.nc
ln -s ${ERAdir}/SST_2020_01_02_ERA5_hourly_Barbados_p021K.nc SST_2020_01_02_ERA5_hourly_Barbados.nc

# SST lower
cd ${currentdir}/../input/perturbations/ERA_SST_m042K
rm SST_2020_01_02_ERA5_hourly_Barbados.nc
ln -s ${ERAdir}/SST_2020_01_02_ERA5_hourly_Barbados_m042K.nc SST_2020_01_02_ERA5_hourly_Barbados.nc

# SST upper
cd ${currentdir}/../input/perturbations/ERA_SST_p042K
rm SST_2020_01_02_ERA5_hourly_Barbados.nc
ln -s ${ERAdir}/SST_2020_01_02_ERA5_hourly_Barbados_p042K.nc SST_2020_01_02_ERA5_hourly_Barbados.nc

# ERA q lower
cd ${currentdir}/../input/perturbations/ERA_q_m30pct
rm q_2020_01_02_ERA5_hourly_Barbados.nc
ln -s ${ERAdir}/q_2020_01_02_ERA5_hourly_Barbados_m30pct.nc q_2020_01_02_ERA5_hourly_Barbados.nc

# ERA q upper
cd ${currentdir}/../input/perturbations/ERA_q_p30pct
rm q_2020_01_02_ERA5_hourly_Barbados.nc
ln -s ${ERAdir}/q_2020_01_02_ERA5_hourly_Barbados_p30pct.nc q_2020_01_02_ERA5_hourly_Barbados.nc

# sonde q lower
cd ${currentdir}/../input/perturbations/sondes_q_m3pct
rm all_sondes.nc
ln -s ${sondedir}/all_sondes_m3pct.nc all_sondes.nc

# sonde q upper
cd ${currentdir}/../input/perturbations/sondes_q_p3pct
rm all_sondes.nc
ln -s ${sondedir}/all_sondes_p3pct.nc all_sondes.nc

echo 'directory structure created'
echo "call ls -l ${currentdir}/../input/perturbations/* to see the input files chosen for each case" 

exit 0
