ifeq ($(USER),bfildier)
	RRTMGP_ROOT = ../rte-rrtmgp
	NFHOME = /usr/local
	NCHOME = /usr/local
endif
ifeq ($(USER),ludo)
	RRTMGP_ROOT = ../../rte-rrtmgp
	NFHOME = /usr/local
	NCHOME = /usr/local
endif

# Location of RTE+RRTMGP libraries, module files.
RRTMGP_BUILD = $(RRTMGP_ROOT)/build
# Sets macros FC, FCFLAGS consistent with RTE+RRTMGP
-include $(RRTMGP_BUILD)/Makefile.conf

# Location of netcdf C and Fortran libraries. Could specify with environment variables if file doesn't exist
-include $(RRTMGP_ROOT)/examples/rfmip-clear-sky/Makefile.libs
#
# RRTMGP library, module files
#
LDFLAGS   += -L$(RRTMGP_BUILD)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(RRTMGP_BUILD)

#
# netcdf library, module files
# C and Fortran interfaces respectively
#
FCINCLUDE += -I$(NFHOME)/include
LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib
LIBS      += -lnetcdff -lnetcdf

VPATH = $(RRTMGP_ROOT)/examples

# Compilation rules
%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)


#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS  = mo_load_coefficients.o mo_simple_netcdf.o

#
# Targets
#
all: sonde_radiation

sonde_radiation: $(ADDITIONS) sonde_radiation.o

sonde_radiation.o: $(ADDITIONS) sonde_radiation.F90

mo_load_coefficients.o: mo_simple_netcdf.o mo_load_coefficients.F90

clean:
	-rm sonde_radiation *.o *.optrpt ../*.optrpt *.mod
