# EUREC4A Radiative profiles

##Steps to run the script

### 0. Add output/ in .gitignore

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

In script, edit compute_radiation_from_soundings.sh with the path to your dropsonde file.
Run compute_radiation_from_soundings.sh

