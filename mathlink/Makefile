
###############################

CC = gcc
CFLAGS = -lgfortran -fopenmp -Wall 

ifdef INTEL
## intel
CC = icc 
CFLAGS = -openmp -Wall -lifcore
ISUFF = _i
endif

###############################

# Directory where the *.o files and the trgfast library are located 
BUILD = ../build$(ISUFF)/
LIBTRG = libtrgfast-0.1.a

#  all: main_c main_f90 trg_link 
all: trg_link 

clean:
	rm -f trg_link main_mma.o main_mmatm.o main_mmatm.c

.PHONY: clean

###############################
# Mathematica related stuff. Adjust paths accordingly.
MLINKDIR = /opt/Wolfram/Mathematica/9.0.1/SystemFiles/Links/MathLink/DeveloperKit
SYS = Linux-x86-64# Set this value with the result of evaluating $SystemID
CADDSDIR = ${MLINKDIR}/${SYS}/CompilerAdditions

INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}

EXTRALIBS = -lm -lpthread -lrt -lstdc++ # Set these with appropriate libs for your system.
MLLIB = ML64i3  # Set this to ML32i3 if using a 64-bit system

MPREP = ${CADDSDIR}/mprep
###############################

# Build the Mathematica things
trg_link: main_mmatm.o main_mma.o $(BUILD)$(LIBTRG) Makefile
	${CC} $(CFLAGS) -o $@ -I${INCDIR} main_mmatm.o main_mma.o -L${LIBDIR} -l${MLLIB} ${EXTRALIBS} $(BUILD)$(LIBTRG)

main_mma.o: main_mma.c 
	$(CC) $(CFLAGS) -c -I$(INCDIR) $< -o $@

main_mmatm.o: main_mmatm.c
	$(CC) $(CFLAGS) -c -I$(INCDIR) $< -o $@
	
main_mmatm.c: main_mma.tm
	${MPREP} $? -o $@
