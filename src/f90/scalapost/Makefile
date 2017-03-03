TARGETS := scalapost

all : $(TARGETS)

scalapost  : scalapost.o 
	$(MPF90) -o $@ $^  ../include/libquietscala.a $(SCALAPACK_LIBS) $(LINK)

scalapost.o : scalapost.f90 ../include/libquietscala.a ../include/libquiet.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ $(TARGETS) 
