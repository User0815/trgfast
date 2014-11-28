# TRG Makefile 
# run `make INTEL=1` for using the Intel compiler suite 

###############################

## gfortran (needs version >=4.6)
F90C = gfortran
# For development:
#  F90FLAGS = -I. -J$(BUILD) -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 \
			-O2 -cpp -fopenmp -g -fcheck=all -fbacktrace
# For production:
#
F90FLAGS = -I. -J$(BUILD) -fPIC  -fmax-errors=1 \
			-O3 -cpp -fopenmp -ffast-math -funroll-loops -march=native 
CC = gcc
CFLAGS = -lm -lgfortran -fopenmp -Wall 

ifdef INTEL
## intel
F90C = ifort 
# For development:
#  F90FLAGS = -module $(BUILD) -fpp2 -O2 -fPIC -ip -CB -vec_report0 -openmp -traceback -g -assume realloc-lhs
# For production:
F90FLAGS = -module $(BUILD) -fpp2 -O3 -CB -openmp -fPIC \
			-assume realloc-lhs -funroll-loops -march=native 
CC = icc
CFLAGS = -lifcore -openmp -Wall 
ISUFF = _i
endif

###############################



##########################
# Do not modify anything below this line.

MODLIST = options tools nonlinear background ode
BUILD_BASE = build
BUILD = $(BUILD_BASE)$(ISUFF)/
SRC = src/
MATHLINK = mathlink/
LIBTRG = libtrgfast-0.1.so
STATLIBTRG = libtrgfast-0.1.a
# all the Fortran IV modules in src/*.f:
FIVMODULES = $(patsubst $(SRC)%.f,%.o,$(wildcard $(SRC)*.f))

MODULES = $(FIVMODULES) $(addsuffix .o, $(MODLIST))

OFILES = $(addprefix $(BUILD),$(MODULES))
  
all: $(BUILD)$(LIBTRG) $(BUILD)$(STATLIBTRG) mathlink driver driver_c

mathlink:
	+$(MAKE) -C $(MATHLINK) 

$(BUILD)%.o: $(SRC)%.f Makefile
	$(F90C) -fPIC -c $< -o $@
	
$(BUILD)%.o: $(SRC)%.f90 Makefile
	$(F90C) $(F90FLAGS) -c $< -o $@
  
$(BUILD)$(STATLIBTRG): $(OFILES) Makefile
	ar rcf $@ $(OFILES) 

$(BUILD)$(LIBTRG): $(OFILES) Makefile
	$(F90C) $(F90FLAGS) -shared -o $@ $(OFILES)

.PHONY: clean mathlink
clean:
	rm -f $(BUILD_BASE)_i/* $(BUILD_BASE)/* driver*
	+$(MAKE) -C $(MATHLINK) clean

driver: $(OFILES) $(SRC)driver.f90
	$(F90C) $(F90FLAGS) $(SRC)driver.f90 -o $@$(ISUFF) $(BUILD)$(STATLIBTRG) 

driver_c: $(OFILES)
	#	$(CC) $(CFLAGS) $(SRC)driver.c -o $@$(ISUFF) $(BUILD)$(STATLIBTRG) 
	+$(MAKE) -C c_wrapper
	cp c_wrapper/driver_c .
