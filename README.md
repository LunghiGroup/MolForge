# MolForge
A multipurpose Fortran2003+MPI suite of codes 
for the simulation of spin dynamics at the 
atomistic level. The code is developed to readily 
interface with first-principles codes.

## Licence and Copyright


## How to Cite MolForge

When using any component of MolForge, please always reference:

- How do phonons relax molecular spins? Science Advances, 5 (9), eaax7163, 2019

- The limit of spin lifetime in solid-state electronic spins, The Journal of Physical Chemistry Letters, 11 (15), 6273-6278, 2020 

- Multiple spin-phonon relaxation pathways in a Kramer single-ion magnet, The Journal of Chemical Physics 153 (17), 174113, 2020

## How to Build MolForge

Molforge is known to compile on Linux machines with GNU 
compilers and openmpi 

MolForge depends on the linear algebra routines Lapack, Blas, and
Scalapack.

To compile and link MolForge and its programs, it 
is necessary to install the command make. 
The available make instances are: all, libMolForge, Spiral, Phondy, 
Tools, clean, clean_MolForge, clean_Spiral, clean_PhonDy,
 clean_Tools.

To compile all codes in one go:

> cd MolForge/

> make all

After successful installation, the folders libs/ and 
bin/ should be present in the molforge main directory. 
The library libMolForge will be build statically and 
move to libs/. All the executables will be move to bin/.
