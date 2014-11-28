CC = gcc
CFLAGS = -lm -lgfortran -fopenmp -Wall -pedantic -g -std=gnu99 -O3

trgfastlib=../build/libtrgfast-0.1.a

all: driver_c

driver_c: driver.c trgfast.o $(trgfastlib)
	$(CC) $(CFLAGS) trgfast.o driver.c -o $@ $(trgfastlib)

trgfast.o: trgfast.c trgfast.h
	$(CC) $(CFLAGS) -c trgfast.c -o $@


.PHONY: clean

clean:
	rm driver_c trgfast.o
