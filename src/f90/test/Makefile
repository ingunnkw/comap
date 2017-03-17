OBJS = time_test.o

all : time_test

time_test : $(OBJS)
	$(MPF90) -o time_test $(OBJS) $(LINK)

time_test.o : time_test.f90 ../include/libquiet.a

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
	@rm -f *.o *.mod *.MOD *.a *~ time_test

