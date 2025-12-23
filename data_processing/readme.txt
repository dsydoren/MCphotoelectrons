This directory contain post-processing and auxiliary programs required to obtain certain data used 
in the manuscript. All programs are in fortran. If the computer has linux and mpif90, one can run
the make command (a Makefile is included). If the fortran compiler is different, one has to edit
the Makefile accordingly or compile the programs separately.

The data processing programs and their description are listed below.

getfluxesfromevdf.out

The source files are igrf13.for and dataproc_get_photoelectron_fluxes_from_evdf_igrf.f90. Note that 
this program includes the IGRF13 model. The IGRF license file from github IGRF_LICENSE.txt is in 
this directory as well. This program uses electron velocity distribution functions calculated by the 
photoelectron transport code to calculate (a) electron distribution functions over energy and pitch
angle and (b) average differential photoelectron fluxes in selected energy ranges. The pitch angle
range for averaging corresponds to the field of view of the photoelectron instrument in the location
of the spacecraft. The numerical grid for the distribution function over energy and pitch angle and
the energy ranges are specified in file input_photoelectron_fluxes_igrf.dat. Another required input
file is input_spacecraft_orbit.dat which provides spacecraft positions at given times. The IGRF and
the orbit data are required to calculate the angle between the axis of the photoelectron instrument
and the geomagnetic field.

getmeasuredspectra.out

The source file is dataproc_get_photoe_fluxes_from_measured_spectra.f90. This program uses measured 
photoelectron spectra to calculate average differential photoelectron fluxes in spacecraft location
versus time. The fluxes correspond to the same energy ranges as defined by the input file
input_photoelectron_fluxes_igrf.dat for the aforementioned program getfluxesfromevdf.out. The output 
of this program is in directory measurements_AEE_orbit_458.

getanglelimitfiles.out

The source file is get_angle_limit_files.f90. This program prepares boxes in the energy-pitch angle 
space and the parallel velocity-transverse velocity space used in Figures 3 and 4 to visualize
limits defining specific photoelectron fluxes. The output of this program is in directory
flux_limits_datafiles.
