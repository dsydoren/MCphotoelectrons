This package contains the input data and software used to prepare manuscript 
"Monte Carlo photoelectron code and its validation using Atmospheric Explorer E data" 
by D. Sydorenko, R. Rankin, J. Liang, and E. Donovan submitted to journal Earth and 
Space Science.

The following directories are included in the package:

calculation_AEE_orbit_458_all_input_files

data_processing

MC_photoelectron_transport_code

Below is the brief decsription of the contents of these directories. Note that each 
directory has a readme.txt file with detailed explanations of all included files.

calculation_AEE_orbit_458_all_input_files -----------------------------------------------

This directory contains all input files necessary to perform calculation of photoelectron 
velocity distribution functions along orbit 458 of AE-E. The input files include both the 
files with simulation parameters that can be modified by the user and the data files required 
by the IRI and MSIS models. This directory can be used to repeat the calculations 
described in the manuscript or to perform simlar calculations for different spacecraft 
locations and times.

data_processing -------------------------------------------------------------------------

This directory contains all post-processing programs.

MC_photoelectron_transport_code ---------------------------------------------------------

This directory contains the source files (fortran and C) of the Monte-Carlo photoelectron 
transport code. The source files include the kinetic photoelectron transport module 
developed by Dmytro Sydorenko at the University of Alberta and the IGRF-13, NRLMSIS 2.1, 
and IRI-2020 models.
