FC := mpif90
CC := mpicc

FFLAGS := -O2
CFLAGS := -O2

OBJ := msis_constants.o \
       msis_utils.o \
       msis_init.o \
       msis_gfn.o \
       msis_tfn.o \
       msis_dfn.o \
       msis_calc.o \
       msis_gtd8d.o \
          cira.o \
          igrf.o \
          iridreg.o \
          iriflip.o \
          irifun.o \
          irirtam.o \
          irisub.o \
          iritec.o \
          rocdrift.o \
       igrf13.o \
          WELL19937a_new.o \
       RandomNumberInterface.o \
          mc_Modules.o \
          mc_PhotoIndividualCollFreqs.o \
          mc_CoordTransform.o \
          mc_EUVInitiateArrays.o \
          mc_EUVproc.o \
          mc_GeneralServiceProc.o \
          mc_MainProgram.o \
          mc_NeutralsMSIS.o \
          mc_PhotoCollFreqs.o \
          mc_PhotoCalculateRates.o \
          mc_PhotoPrepareFieldline.o \
          mc_PhotoProc.o

runphotoe : $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o : IRI-2020-source/%.for
	$(FC) $(FFLAGS) -c $< 

%.o : nrlmsis2.1/%.F90 
	$(FC) $(FFLAGS) -c $< 

%.o : MC_transport/%.for
	$(FC) $(FFLAGS) -c $< 

%.o : MC_transport/%.f90
	$(FC) $(FFLAGS) -c $< 

%.o : MC_transport/%.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.mod *.o runphotoe
