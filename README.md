# MolForge
A multipurpose FORTRAN2003 and MPI suite of codes for the simulation of quantum spin dynamics at the atomistic level. The code is developed to readily 
interface with electronic structure codes and makes it possible to predict spin decoherence and spin-phonon relaxation fully ab initio.

## Licence and Copyright

MolForge: a multipurpose FORTRAN2003 and MPI suite of codes for the simulation of quantum spin dynamics at the atomistic level. 

Copyright (c) 2021-2024 Trinity College Dublin

Developed by Alessandro Lunghi, School of Physics, AMBER and CRANN Institute, Trinity College, Dublin 2, Ireland

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

## How to Cite MolForge

When using  MolForge or any of its components, please always cite: 

- Toward exact predictions of spin-phonon relaxation times: An ab initio implementation of open quantum system theory, Science Advances 8(31), eabn7889, 2022

Additional works that contributed towards the formalism implemented in MolForge are:

- How do phonons relax molecular spins? Science Advances, 5 (9), eaax7163, 2019

- Electronic spin-spin decoherence contribution in molecular qubits by quantum unitary spin dynamics, Journal of Magnetism and Magnetic Materials 487, 165325, 2019

- The limit of spin lifetime in solid-state electronic spins, The Journal of Physical Chemistry Letters, 11 (15), 6273-6278, 2020 

- Multiple spin-phonon relaxation pathways in a Kramer single-ion magnet, The Journal of Chemical Physics 153 (17), 174113, 2020

## How to Build MolForge

Molforge is known to compile on Linux machines with GNU compilers(v12.2.0) and openmpi(5.0.2).

MolForge depends on the linear algebra routines Lapack (v3.10.0), and Netlib Scalapack(2.2.0).

The following variables in the main Makefile might need to be set 

> export INC_DIR = "path to libraries include files"

> export LIB_DIR = "path to libraries files"

> export LIBS = -lscalapack -llapack -lblas

> export INC_FLAG = -I$(INC_DIR) 
 
> export LINK_FLAG = -L$(LIB_DIR) 


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
