l3gen : l3gen.o ../include/libquiet.a
	$(MPF90) -o $@ $< $(LINK)

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	icc $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ l3gen

