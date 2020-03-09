# This disables builtin rules. Important to avoid surprising breakage with
# pattern rules
.SUFFIXES:

OBJS  := quiet_defs.o quiet_healpix_mod.o quiet_target_mod.o \
         quiet_pixspace_mod.o quiet_module_mod.o comap_detector_mod.o \
         comap_frequency_mod.o comap_ephem_mod.o comap_acceptlist_mod.o comap_scan_mod.o \
         comap_jackknife_mod.o d1mach.o drc3jj.o quiet_l1_defs.o \
         math_tools.o powell_mod.o quasi_newton_mod.o cholesky_decomposition_mod.o mjd_to_gregorian_mod.o \
         quiet_fileutils.o quiet_detector_mod.o quiet_mapfile_mod.o \
         quiet_nr_mod.o quiet_utils.o quiet_hdf_mod.o quiet_Lx_mod.o \
	 comap_defs.o comap_lx_mod.o comap_ces_mod.o comap_map_mod.o \
         comap_gain_mod.o spline_1D_mod.o spline_2D_mod.o sort_utils.o \
	 quiet_mpi_mod.o fake_mpi_mod.o quiet_pointing_mod.o comap_pointing_mod.o \
         quiet_apex_mod.o quiet_gain_mod.o quiet_mpi_utils.o \
         quiet_sidelobe_mod.o quiet_fft_mod.o \
         quiet_i2qu_mod.o quiet_mask_mod.o quiet_pmac_mod.o \
         quiet_assembly_mod.o quiet_filter_mod.o quiet_postutils.o \
         quiet_noise_estimation_mod.o comap_patch_mod.o quiet_todsim_mod.o \
         quiet_glitchy_noise_mod.o quiet_task_mod.o comap_noise_estimation_mod.o \
         quiet_system_mod.o quiet_system_backend.o \
         quiet_shared_output_mod.o ziggurat.o sdeOOF.o quiet_status_mod.o \
         quiet_stat_mod.o quiet_ephem_backend.o locate_mod.o \
         quiet_constrained_mod.o comap_patch_detect_mod.o quiet_typeb_mod.o \
         quiet_pointing_mod2.o quiet_covfile_mod.o quiet_angcorr_mod.o  \
         quiet_ephem_mod.o ../../cpp/ephcom/ephcom.o ../../cpp/ephcom/gnulliver.o wmap_beam_mod.o \

SCALAOBJS := scalawrap.o quiet_mpi_mod.o scalautils.o

libquiet.a : $(OBJS)
	$(AR) $(ARFLAGS) libquiet.a $(OBJS)
	$(RANLIB) libquiet.a

libquietscala.a: $(SCALAOBJS)
	$(AR) $(ARFLAGS) $@ $(SCALAOBJS)
	$(RANLIB) $@

quiet_angcorr_mod.o      : spline_1D_mod.o
quiet_covfile_mod.o      : quiet_hdf_mod.o
quiet_pointing_mod2.o    : quiet_detector_mod.o quiet_utils.o
quiet_mpi_utils.o        : quiet_mpi_mod.o quiet_utils.o quiet_system_mod.o
quiet_pmac_mod.o         : quiet_utils.o locate_mod.o quiet_fileutils.o
quiet_detector_mod.o     : quiet_utils.o
comap_detector_mod.o     : quiet_utils.o
comap_frequency_mod.o    : quiet_utils.o
quiet_mask_mod.o         : quiet_utils.o quiet_fileutils.o
comap_patch_detect_mod.o : comap_patch_mod.o comap_detector_mod.o quiet_pointing_mod.o
comap_patch_mod.o        : quiet_utils.o quiet_ephem_mod.o
quiet_constrained_mod.o  : quiet_utils.o quiet_fft_mod.o math_tools.o
quiet_ephem_mod.o        : comap_pointing_mod.o
quiet_healpix_mod.o      : quiet_utils.o quiet_system_mod.o
quiet_pixspace_mod.o     : quiet_utils.o quiet_fileutils.o math_tools.o spline_1D_mod.o
quiet_shared_output_mod.o : quiet_system_mod.o quiet_utils.o quiet_mpi_mod.o
quiet_status_mod.o       : quiet_system_mod.o quiet_utils.o quiet_shared_output_mod.o
quiet_task_mod.o         : quiet_system_mod.o
quiet_sidelobe_mod.o     : quiet_utils.o quiet_fileutils.o quiet_module_mod.o
quiet_system_mod.o       : quiet_system_backend.o
quiet_acceptlist_mod.o   : quiet_ces_mod.o quiet_module_mod.o
comap_acceptlist_mod.o   : comap_scan_mod.o comap_detector_mod.o comap_frequency_mod.o
quiet_target_mod.o       : quiet_utils.o quiet_ces_mod.o quiet_acceptlist_mod.o quiet_defs.o
#l1_read_comap_mod.o      : quiet_utils.o 
comap_scan_mod.o         : quiet_utils.o
comap_ephem_mod.o        : spline_1D_mod.o
quiet_ces_mod.o          : quiet_utils.o
quiet_Lx_mod.o           : quiet_hdf_mod.o quiet_fft_mod.o
comap_lx_mod.o		 : quiet_mpi_mod.o quiet_hdf_mod.o quiet_fft_mod.o
comap_map_mod.o          : quiet_utils.o quiet_hdf_mod.o
quiet_filter_mod.o       : quiet_noise_estimation_mod.o quiet_assembly_mod.o
quiet_fileutils.o        : quiet_utils.o quiet_hdf_mod.o quiet_mapfile_mod.o quiet_covfile_mod.o
quiet_nr_mod.o           : math_tools.o
spline_1D_mod.o          : math_tools.o locate_mod.o
spline_2D_mod.o          : spline_1D_mod.o
quiet_utils.o            : sort_utils.o spline_1D_mod.o quiet_system_mod.o quiet_defs.o
quiet_i2qu.o             : quiet_fileutils.o pix_tools.o
quiet_assembly_mod.o     : quiet_utils.o
scalawrap.o              : quiet_mpi_mod.o
quiet_pointing_mod.o     : quiet_utils.o math_tools.o quiet_pmac_mod.o
comap_pointing_mod.o     : quiet_utils.o math_tools.o 
scalautils.o             : quiet_mpi_mod.o scalawrap.o
quiet_noise_estimation_mod.o : quiet_utils.o quiet_fft_mod.o 
comap_noise_estimation_mod.o : quiet_utils.o quiet_fft_mod.o 
quiet_module_mod.o       : quiet_utils.o quiet_detector_mod.o
quiet_todsim_mod.o       : quiet_gain_mod.o quiet_assembly_mod.o ziggurat.o
quiet_stat_mod.o         : math_tools.o
quiet_glitchy_noise_mod.o : quiet_utils.o quiet_constrained_mod.o
wmap_beam_mod.o           : quiet_utils.o
comap_jackknufe_mod.d     : quiet_utils.o quiet_hdf_mod.o

quiet_system_backend.o: quiet_system_backend.cpp
	$(MPCC) $(CCOMP) -DREVISION=$$(svn info | awk '/^Revis[ij]on/{print $$2}') -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) $(OPENMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

%.o : %.cpp
	$(MPCC) $(CCOMP) -I../../cpp/ephcom -c $< 

%.f90 : %.f90.in
	$(TEMPITA) < "$<" > "$@"

clean :
	@rm -f *.o *.mod *.MOD *.a *~

