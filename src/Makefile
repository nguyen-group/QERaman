# Makefile for phonon (PH)

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../dft-d3 \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a
QEMODS = $(BASEMODS)
TLDEPS= phlibs pwlibs
PHOBJS = phq_readin2.o do_phonon2.o supp.o ../../PHonon/PH/libph.a
PPOBJS = write_p_avg.o ../../PP/src/libpp.a
MODULES = $(PWOBJS) $(QEMODS)

RAOBJS = io_raman.o raman_mod.o knum.o raman_read.o raman_output.o calc_raman_intensity.o
#RAMODS = io_raman.mod raman_mod.mod raman_read.mod raman_output.mod calc_raman_intensity.mod

all : ph_mat.x bands_mat.x raman.x
#all : raman.x

ph_mat.x : phonon2.o $(PHOBJS) $(PWOBJS) $(LRMODS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ phonon2.o $(PHOBJS) \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( mkdir ../bin; cd ../bin ; ln -fs ../src/ph_mat.x . )

bands_mat.x : bands.o $(PPOBJS) $(MODULES)
	$(LD) $(LDFLAGS) -o $@ bands.o $(PPOBJS) \
		$(MODULES) $(QELIBS)
	- ( cd ../bin ; ln -fs ../src/bands_mat.x . )

raman.x : raman.o $(RAOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ raman.o $(RAOBJS) \
		$(QEMODS)
	- ( cd ../bin ; ln -fs ../src/raman.x . )

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L ../bin/*

depend:
	@echo 'Checking dependencies...'
	- ( if test -x install/makedeps.sh ; then install/makedeps.sh ; fi)

phonon2.o: supp.o
do_phonon2.o: supp.o
phq_readin2.o: supp.o
bands.o: write_p_avg.o
raman.o: io_raman.o
raman.o: raman_mod.o
raman.o: knum.o
raman.o: raman_read.o
raman.o: raman_output.o
raman.o: calc_raman_intensity.o

include ../../PP/src/make.depend ../../PHonon/PH/make.depend
