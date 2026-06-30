This directory contains all input files required to perform simulation of orbit 458 of AE-E 
and do the post-processing. These include input files with parameters that can be specified 
by the user (input*dat) and files required by the IRI and MSIS models. All user-defined 
input files include text comments describing the input parameters.

--------------- User-defined simulation input files: -----------------------------------------

File input_photoelectrons.dat specifies major simulation parameters.

File input_spacecraft_orbit.dat contains times and coordinates of the 15 orbit points.

File input_photoelectron_fluxes_igrf.dat is used by data processing program getfluxesfromevdf.out,
see readme.txt in directory data_processing.

The following input files enable detailed model diagnostics. The calculation can proceed without
these files. Note that enabling the diagnostics slows down the calculation. For faster calculation, the
input files mentioned below can be removed from the working directory.

input_coll_freq_ranges_op_PPP.dat

These files specify energy ranges for calculation of average electron-neutral collision frequencies.
In the file name, the PPP is the 3-digit number of the orbit point (001, 002, ..., 015). If the program
finds such a file, it calculates frequencies of all electron-neutral collisions averaged over the
requested energy ranges along the geomagnetic field line through the requested orbit point. Theoretical
values of frequencies are saved in files op_PPP_collfreq_s1_avgw_EEE_eV_pfl.dat, the frequencies
actually measured in the simulation are saved in files op_PPP_collfreq_s1_avgw_EEE_eV_measured_pfl.dat.
Here PPP is the orbit point number, EEE is the middle energy of the averaging energy range. These
output files have text headers explaining the content.

Files with suffix "pfl" save data along the whole field line. There are also files containing the
corresponding data for the spacecraft location only, with names op_PPP_collfreq_s1_avgw_EEE_eV_osp.dat
and op_PPP_collfreq_s1_avgw_EEE_eV_measured_osp.dat, respectively (note that "pfl" is replaced with
"osp"). Suffixes "pfl" and "osp" stand for "photoelectrons [along the] field line" and
"only spacecraft [location]", respectively

input_energy_pa_ranges_op_PPP.dat

These files specify energy and pitch angle limits for photoelectron fluxes in the selected orbit points.
In the file name, the PPP is the orbit point number. If the program finds such a file, it calculates
the photoelectron fluxes in the requested energy ranges and the pitch angle range corresponding to the
field of view of the AE-E photoelectron instrument. The fluxes are calculated along the geomagnetic
field line through the requested orbit point. The program also calculates all possible sources and
sinks contributing to these fluxes. The output is saved in files op_PPP_range_RR_pfl.dat, where the PPP
is the orbit point number and RR is the two-digit energy range number. Files op_PPP_range_RR_pfl.dat
have a text header explaining the file content. The program also saves the values of the same dataset
as in op_PPP_range_RR_pfl.dat in the spacecraft location only, these values are saved in files
op_PPP_range_RR_osp.dat.


