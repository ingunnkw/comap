TARGETS = fsmap cov_convert rect_convolve season_weight template_fit gainfit2 corrstat gainstat gainfit chisqstat jackstat pixaccept freq2corr pixtrans hunf fullval qmount tilt2hdf tilt2ces coord_convert pixsmooth addpix_cov addrhs_eqn region_split region_paint covslice pixdeg fits2cov l1_mjd euler sunmask lxcat lxmulticat season_corr convert_covar_mat eqn_convert fitsimg wtemplate #fixtone

default: $(TARGETS)

% : %.o ../include/libquiet.a
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
	@rm -f *.o *.mod *.MOD *.a *~ $(TARGETS)

