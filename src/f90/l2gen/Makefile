all: l2gen

l2gen : l2gen.o ../include/libquiet.a
	$(MPF90) -o $@ $< $(LINK)

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

clean :
	@rm -f *.o *.mod *.MOD *.a *~ l2gen

