Requirements
============

You need a Fortran compiler supporting the Fortran 2008 standard.
Successfully tested were the compiler included in the gcc v4.6.2, and the
Intel Compiler Suite 12.0.2.  

When compiling with `gcc`, make sure both `gfortran` and `gcc` are >= v4.6.2
and that the appropriate libraries are in your `LD_LIBRARY_PATH` variable.

If you plan to call the library from within Wolfram Mathematica via
MathLink, you obviously need Mathematica. You may have to edit the Makefile
in the `mathlink` directory to set the paths according to your system.



How to use the code
===================

Type `make` to build the library with the `gcc` compilers, and `make
INTEL=1` to build it with the Intel compilers. It can then be called from
your Fortran or C/C++ code. You'll find examples on how to use it in the
files `src/driver.f90` and `src/driver.c`. A demo notebook that shows how to
use the MathLink is located at `mathlink/trgfast-demo.nb`.

To build your own program using the trgfast library, use this command to
compile when using Fortran

	gfortran -fopenmp my_program.f90 libtrgfast-0.1.a

and this when using C

    gcc -lgfortran -fopenmp my_program.c libtrgfast-0.1.a

You may have to adjust the path to the file `libtrgfast-0.1.a` or copy it
into the same directory.

Troubleshooting
===============

* The MathLink died? If there are no messages on the command line, try
  making sure that the trgfast library and the `trg_link` binary are compiled
  on the same physical machine as the one that is running Mathematica.

Feedback
========

Please send bug reports, complaints, and other feedback to
vollmer@thphys.uni-heidelberg.de.
