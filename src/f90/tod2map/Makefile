all : tod2comap covbench comaptest #maptest

tod2comap : libtod2map.a tod2comap.o ../include/libquiet.a
	$(MPF90) -o tod2comap tod2comap.o libtod2map.a $(LINK) $(COMMON_FLAGS)
libtod2map.a : tod2comap_utils.o tod2comap_cl_mod.o tod2comap_mapmaker.o
	ar crs $@ $^
covbench : libtod2map.a covbench.o ../include/libquiet.a
	$(MPF90) -o covbench covbench.o libtod2map.a $(LINK) $(COMMON_FLAGS)
maptest : libtod2map.a maptest.o ../include/libquiet.a
	$(MPF90) -o $@ maptest.o libtod2map.a $(LINK) $(COMMON_FLAGS)
comaptest : comaptest.o 
	$(MPF90) -o comaptest comaptest.o $(LINK) $(COMMON_FLAGS)

tod2comap_cl_mod.o   : tod2comap_utils.o
tod2comap_mapmaker.o : tod2comap_utils.o

%.o : %.f90
	$(MPF90) $(F90COMP) $(OPENMP) -c $< -Os

clean :
	@rm -f *.o *.mod *.MOD *.a *~ tod2comap

