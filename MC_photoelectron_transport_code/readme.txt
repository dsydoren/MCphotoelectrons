This software calculates distribution functions of photoelectrons over the parallel and
transverse velocities (relative to the geomagnetic field) in given locations and times. It is
expected, though it is not critical, that the locations and times correspond to points along
an orbit of a spacecraft. The model considers closed field lines if the apex altitude of the
field line does not exceed 2000 km (due to the upper limit of altitudes considered by the IRI
model). Otherwise, the model treats the field line as an open one.

The software includes the following:

The photoelectron transport module developed at the University of Alberta by Dmytro Sydorenko
consisting of files mc_*f90 in directory MC_transport.

The NRLMSIS 2.1 model of neutral atmosphere in directory nrlmsis2.1. The readme.txt and the
license file nrlmsis2.1_license.txt included with the NRLMSIS 2.1 software are in directory
nrlmsis2.1.

The International Geomagnetic Reference Field IGRF-13 in file igrf13.for in directory
MC_transport. The IGRF license file IGRF_LICENSE.txt downloaded from
https://github.com/space-physics/igrf/blob/main/LICENSE.txt is in directory MC_transport.

The International Reference Ionosphere Software (Dec 5, 2023) in directory IRI-2020-source.
The 00readme.txt and the license file 00_iri-License.txt downloaded from
https://irimodel.org/IRI-2020/ are in directory IRI-2020-source.

The random number generator by Shin Harase, Hiroshima University, which is a "maximally
equidistributed" version of the random number generator WELL19937a by Francois Panneton and
Pierre L'Ecuyer, University of Montreal and Makoto Matsumoto, Hiroshima University, in file
WELL19937a_new.c in directory MC_transport.

The fortran interface for the random number generator developed by Salomon Janhunen (Los
Alamos National Laboratory) in file RandomNumberInterface.f90 in directory MC_transport.

The data files required by the IRI and MSIS models, as well as the input files with parameters
that can be adjusted by the user are in directory
../calculation_AEE_orbit_458_all_input_files.

The software requires a fortran compiler with an MPI library and a C compiler. The Makefile
included with the software assumes that the MPI fortran compiler is mpif90 and the C compiler
is mpicc. If your computer has these compilers, it is sufficient to run command make in the
same directory as the Makefile. Otherwise one has to edit the Makefile accordingly.

The software is supposed to run on multiple CPU cores. To run the software, copy the compiled
executable file (default name runphotoe) into the directory with the input files and run the
following command (for example, if the computer has 16 CPU cores):

mpirun -np 16 ./runphotoe

or

mpirun -np 16 ./runphotoe > output.txt &

or 

nohup mpirun -np 16 ./runphotoe &

The user-adjusted input files and the output files are discussed in readme.txt in directory
../calculation_AEE_orbit_458.





