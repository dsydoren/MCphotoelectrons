This directory contains all input files required to perform simulation of orbit 458 of AE-E 
and do the post-processing. These include input files with parameters that can be specified 
by the user (input*dat) and files required by the IRI and MSIS models. All user-defined 
input files include text comments describing the input parameters.

--------------- User-defined simulation input files: -----------------------------------------

File input_photoelectrons.dat specifies major simulation parameters.

File input_spacecraft_orbit.dat contains times and coordinates of the 15 orbit points.

Files input_energy_pa_ranges_op_*.dat specify energy and pitch angle limits for photoelectron 
fluxes in the 15 orbit points. In the file name, the * is the 3-digit number of the orbit 
point (001, 002, ..., 015). These are diagnostics files which are not necessary to obtain 
the velocity distribution functions. The diagnostics output (saved in files op_*_pfl.dat),
is used in Figures 6, 7, D1, and D2.



