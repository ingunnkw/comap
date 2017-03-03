OBJS = maptool.o

all : maptool

maptool  : $(OBJS)
	$(MPF90) -o maptool $(OBJS) $(LINK)

maptool.o : maptool.f90 ../include/libquiet.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ 

