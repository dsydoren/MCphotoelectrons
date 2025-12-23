FC := mpif90

FFLAGS := -O2

all: getanglelimitfiles.out getfluxesfromevdf.out getmeasuredspectra.out

OBJ0 := get_angle_limit_files.f90

OBJ1 := igrf13.for dataproc_get_photoelectron_fluxes_from_evdf_igrf.f90

OBJ2 := dataproc_get_photoe_fluxes_from_measured_spectra.f90

getanglelimitfiles.out : $(OBJ0)
	$(FC) $(FFLAGS) -o $@ $(OBJ0)

getfluxesfromevdf.out : $(OBJ1)
	$(FC) $(FFLAGS) -o $@ $(OBJ1)

getmeasuredspectra.out : $(OBJ2)
	$(FC) $(FFLAGS) -o $@ $(OBJ2)

clean:
	rm -f getanglelimitfiles.out getfluxesfromevdf.out getmeasuredspectra.out

