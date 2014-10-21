Requirements
============

You need a Fortran compiler supporting the Fortran 2008 standard.
Successfully tested were the compiler included in the gcc v4.6.2, and the
Intel Compiler Suite 12.0.2.  

When compiling with gcc, make sure both `gfortran` and `gcc` are >= v4.6.2
and that the appropriate libraries are in your `LD_LIBRARY_PATH` variable.

If you plan to call the library from within Wolfram Mathematica via
MathLink, you obviously need Mathematica. You may have to edit the Makefile
in the `mathlink` directory to set the paths according to your system.



How to use the code
===================

Type `make` to build the library. It can then be called from your Fortran or
C/C++ code. You'll find examples on how to use it in the files `src/driver.f90`
and `src/driver.c`.



Feedback
========

Please send bug reports, complaints, and other feedback to
vollmer@thphys.uni-heidelberg.de.
