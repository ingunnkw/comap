
COBJS  := map_editor_utils.o map_editor_simple_ops_mod.o map_editor_complex_ops_mod.o 

map_editor : libmap_editor.a map_editor.o
	$(MPF90) -o map_editor map_editor.o libmap_editor.a $(LINK) 

map_editor_complex_ops_mod.o    : map_editor_utils.o 

libmap_editor.a : $(COBJS)
	$(AR) $(ARFLAGS) libmap_editor.a $(COBJS)
	$(RANLIB) libmap_editor.a

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
	@rm -f *.o *.mod *.MOD *.a *~ map_editor

