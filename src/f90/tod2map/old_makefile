all : tod2map covbench #maptest

tod2map : libtod2map.a tod2map.o ../include/libquiet.a
	$(MPF90) -o tod2map tod2map.o libtod2map.a $(LINK) $(COMMON_FLAGS)
libtod2map.a : tod2map_utils.o tod2map_cl_mod.o tod2map_mapmaker.o
	ar crs $@ $^
covbench : libtod2map.a covbench.o ../include/libquiet.a
	$(MPF90) -o covbench covbench.o libtod2map.a $(LINK) $(COMMON_FLAGS)
maptest : libtod2map.a maptest.o ../include/libquiet.a
	$(MPF90) -o $@ maptest.o libtod2map.a $(LINK) $(COMMON_FLAGS)

tod2map_cl_mod.o   : tod2map_utils.o
tod2map_mapmaker.o : tod2map_utils.o

%.o : %.f90
	$(MPF90) $(F90COMP) $(OPENMP) -c $< -Os

clean :
	@rm -f *.o *.mod *.MOD *.a *~ tod2map

