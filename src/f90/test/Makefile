OBJS = #bin_test.o # time_test.o

all : time_test bin_test

time_test : $(OBJS) time_test.o
	$(MPF90) -o time_test $(OBJS) time_test.o $(LINK)

bin_test : $(OBJS) bin_test.o
	$(MPF90) -o bin_test $(OBJS) bin_test.o $(LINK)

time_test.o : time_test.f90 ../include/libquiet.a

bin_test.o : bin_test.f90 ../include/libquiet.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $<

%.o : %.cpp
	$(MPCC) $(CCOMP) -c $<

clean :
	@rm -f *.o *.mod *.MOD *.a *~ time_test bin_test

