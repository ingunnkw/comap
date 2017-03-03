all : tester

tester : tester.o ../include/libquiet.a
	$(MPF90) -o tester tester.o $(LINK) 
#tester : testmod.o tester.o ../include/libquiet.a ../include/libquietscala.a
#	$(MPF90) -o tester tester.o testmod.o ../include/libquietscala.a $(SCALAPACK_LIBS) $(LINK) 

%.o : %.f90
	$(MPF90) $(F90COMP) $(OPENMP) -I../map2cl/ -c $<

clean :
	@rm -f *.o *.mod *.MOD *.a *~ tester

