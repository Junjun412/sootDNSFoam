include Makefile.in

F90FILES = precision.f90 string.f90 soot_mod.f90 soot.f90 soot_chem.f90 soot_momic.f90 soot_hmom.f90 soot_dqmom.f90 soot_cqmom.f90 soot_semi.f90 pda.f90 soot_driver.f90
LIBFILE = libsoot.a

OFILES = $(F90FILES:.f90=.o)
MODFILES = $(F90FILES:.f90=.mod)

.SUFFIXES: .o .f90

$(LIBFILE):$(OFILES) Makefile	
	cd $(OBJDIR); $(AR) $(LIBDIR)/$(LIBFILE) $(OFILES); $(RL) $(LIBDIR)/$(LIBFILE)

.f90.o: 
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) $(LAPACK_LIB) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)

clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(LIBDIR); rm -f $(LIBFILE)
