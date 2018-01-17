COBJS = scan_validate_mod.o

all : scan_validate

scan_validate : libces_validate.a scan_validate.o ../include/libquiet.a
	$(MPF90) -o scan_validate scan_validate.o libces_validate.a $(LINK)

scan_validate.o   : ../include/libquiet.a libces_validate.a

libces_validate.a : $(COBJS)
	$(AR) $(ARFLAGS) libces_validate.a $(COBJS)
	$(RANLIB) libces_validate.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	icc $(CCOMP) -c $<

clean :
	@rm -f *.o *.mod *.MOD *.a *~ ces_validate

