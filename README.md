# M-dwarf Stellar Winds

This repository is where I am developing the Fortran code to calculate
all of the initial conditions for a large-scale simulation. 

It consists of an input file, which contains the three parameters 
needed to constrain the simulations. I.e. the variables to 
adjust in order to get the actuall initial conditions. The 
specifics of this are described in the project file.

The script consists of a few distinct parts; a file
reader that reads the three parameters, the subroutines 
that actually perform the calculations, and the main program that
calls the subroutines.

When run, the program will read the input file, perform all
calculations, and write two output files. The first is a file
that lists the physical parameters of the simulation that will 
be run and the other is the main input file for the large-
scale simulation with all of the input parameters converted
to unitless numbers (which is what the simulation requires).

The first number (integer) that is read in the input file 
is used to keep track of what situations are used. This is 
a good way to make sure that your physical parameter file 
and simulation input match. Whenever you "create" a new
simulation input file, change the first number. That number
will appear in the file names of both the simulation input
file and the physical parameter file, which makes it easy to 
keep up with all of the different runs that will be performed. 

The python version of this script (what was originally used)
will also be included in this repository at some point.

For questions regarding this program, please contact me at...

jamesmccord14@gmail.com



