COBJS  := map2cl_utils.o map2cl_data_mod.o map2cl_nr_mod.o map2cl_parfit_mod.o

all : map2cl

map2cl_data_mod.o   : map2cl_utils.o ../include/libquietscala.a ../include/libquiet.a
map2cl_nr_mod.o     : map2cl_data_mod.o
map2cl_parfit_mod.o : map2cl_nr_mod.o
map2cl.o            : map2cl_parfit_mod.o 
map2cl              : ../include/libquiet.a map2cl_nr_mod.o map2cl.o

map2cl  : libmap2cl.a ../include/libquiet.a ../include/libquietscala.a map2cl.o
	$(MPF90) -o map2cl map2cl.o libmap2cl.a ../include/libquietscala.a $(SCALAPACK_LIBS) $(LINK)

libmap2cl.a : $(COBJS)
	$(AR) $(ARFLAGS) libmap2cl.a $(COBJS)
	$(RANLIB) libmap2cl.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ map2cl

