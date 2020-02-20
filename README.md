# EUREC4A Radiative profiles


## input

Where the profile for a standard tropical atmosphere is stored.

##Steps to run the script

### 0. Add input/Data_local in .gitignore

### 1. Compile the radiation code

In rte-rrtmgp/build create a Makefile.conf following the template in the folder.
Alternatively, you can set environment variables FC and FCFLAGS to be the name of the
Fortran compiler and the compilation flags.
Invoke make.

### 2. Compile main script

In Makefile edit the flags NCHOME and NFHOME for your platform. These point to the
root of the netCDF C and Fortran installations on your platform. 
Call make.

### 3. Run the script

In script, call compute_radiation_from_soundings.sh [full path for the dropsonde netcdf file].
