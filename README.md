# EUREC4A Radiative profiles


## input

Where the profile for a standard tropical atmosphere is stored.

##Steps to run the script


### 1. Compile the radiation code 

In rte-rrtmgp/build create a Makefile.conf following the template in the folder.
Call make.

### 2. Compile main script

In code, edit the flags NCHOME and NFHOME for your platform.
Call make.

### 3. Run the script

In script, call compute_radiation_from_soundings.sh [full path for the dropsonde netcdf file].


