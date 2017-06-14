all: scan_detect

scan_detect : scan_detect.o ../include/libquiet.a
	$(MPF90) -o scan_detect scan_detect.o $(LINK)

#l1_dump : l1_dump.o ../include/libquiet.a
#	$(MPF90) -o $@ $< $(LINK)

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

clean :
	@rm -f *.o *.mod *.MOD *.a *~ scan_detect
