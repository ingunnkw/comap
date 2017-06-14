#
# Adapted from SVN setup written by Ted Kisner
#
# QUIET Makefile for use with EITHER the planck
# module build system OR a set of stand-alone
# config files.  Default make target prints help.
#
# This makefile requires the use of GNU make or
# compatible tools that allow the substitution
# operator ":=".
#
#

NOW := $(shell date +%Y%m%d-%H%M)
DIR := quiet_$(NOW)

TOPDIR := $(shell pwd)
export TOPDIR
export PYTHONPATH := "$(TOPDIR)/src/python:$(PYTHONPATH)"

ifdef TARGET
	ifndef PREFIX
		PREFIX := ./$(TARGET)
	endif
	INSTALL := $(PREFIX)/quiet
else
	ifdef QUIET
		include $(TOPDIR)/config/config.$(QUIET)
		ifndef INSTALL
	INSTALL := $(TOPDIR)/install_$(QUIET)
		endif
	else
		$(error QUIET undefined)UNDEFINED
	endif
	ifndef MAKE
		export MAKE := make
	endif
	ifndef AR
		export AR := ar
	endif
	ifndef ARFLAGS
		export ARFLAGS := crs
	endif
	ifndef RANLIB
		export RANLIB := ranlib
	endif
	ifndef CTAR
		export CTAR := tar czvf
	endif
	ifndef F90
		export F90 := f90
	endif
	ifndef MPF90
		export MPF90 := mpif90
	endif
	ifndef F90FLAGS
		export F90FLAGS := -g -O2
	endif
	ifndef MPF77
		export MPF77 := mpif77
	endif
	ifndef FFLAGS
		export FFLAGS := -g -O2
	endif
	ifndef MPCC
		export MPCC := cc
	endif
	ifndef CFLAGS
		export CFLAGS := -O2
	endif
	ifndef LDFLAGS
		export LDFLAGS := -lm
	endif
	ifndef FORTRAN_UPPER
		export FORTRAN_UPPER := 0
	endif
	ifndef CFITSIO_LINK
		export CFITSIO_LINK := -L/usr/local/lib -lcfitsio
	endif
	ifndef LAPACK_LINK
		export LAPACK_LINK := -L/usr/local/lib -llapack -lblas
	endif
	ifndef F90FLAGS_DEBUG
		export F90FLAGS_DEBUG := -g -O0
	endif
	ifndef FFLAGS_DEBUG
		export FFLAGS_DEBUG := -g -O0
	endif
	ifndef CFLAGS_DEBUG
		export CFLAGS_DEBUG := -g -O0
	endif
	ifndef LDFLAGS_DEBUG
		export LDFLAGS_DEBUG := -lm
	endif
endif

ifdef DEBUG
	export F90FLAGS := $(F90FLAGS_DEBUG)
	export FFFLAGS  := $(FFFLAGS_DEBUG)
	export CFLAGS   := $(CFLAGS_DEBUG)
	export LDFLAGS  := $(LDFLAGS_DEBUG)
endif

export CCOMP := $(CFLAGS)  -I$(TOPDIR)/src/f90/include $(HEALPIX_INCLUDE) $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(FFTW_INCLUDE)

export F90COMP := $(F90FLAGS) -I$(TOPDIR)/src/f90/include $(HEALPIX_INCLUDE) $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(FFTW_INCLUDE) $(HDF_INCLUDE)
export FCOMP := $(FFLAGS) -I$(TOPDIR)/src/f90/include $(LAPACK_INCLUDE) $(HEALPIX_INCLUDE) $(CFITSIO_INCLUDE) $(FFTW_INCLUDE)
export LINK := -L$(TOPDIR)/src/f90/include -lquiet $(HEALPIX_LINK) $(CFITSIO_LINK) $(LAPACK_LINK) $(FFTW_LINK) $(NOVAS_LINK) $(LDFLAGS) $(HDF_LINK) $(LDFLAGS) $(OPENMP)
export TEMPITA := "$(TOPDIR)/src/python/tempita_proc.py"

all : libquiet libutil l2gen postmap map_editor scan_validate maptool test # ces_detect l3gen tod2map utils_f90

full : all libquietscala scalapost map2cl


help :
	@echo ' '
	@echo '  This Makefile is used to build Quiet in a way that is'
	@echo '  customized to your system.  You must export the QUIET'
	@echo '  environment variable and set it to the name of your platform.'
	@echo '  Then you must create a config file named config/config.<$QUIET>.'
	@echo '  I suggest copying the example config file and modifying the'
	@echo '  parameters for your platform.'
	@echo ' '
	@echo '  The following make targets are supported:'
	@echo ' '
	@echo '    make         : build everything'
	@echo '    make help    : print this help screen'
	@echo '    make install : install everything'
	@echo '    make clean   : remove build files'
	@echo '    make dist    : construct a date-stamped source tarball'
	@echo ' '

install : all
	@mkdir -p $(INSTALL)/lib
	@mkdir -p $(INSTALL)/include
	@mkdir -p $(INSTALL)/bin
	@cp src/f90/include/libquiet.a $(INSTALL)/lib
	@cp src/f90/postmap/postmap $(INSTALL)/bin
	@cp src/f90/map2cl/map2cl $(INSTALL)/bin
	@cp src/f90/scalapost/scalapost $(INSTALL)/bin
	@cp src/f90/map_editor/map_editor $(INSTALL)/bin
	@cp src/f90/tod2map/tod2map $(INSTALL)/bin
	@cp src/cpp/utils/{map2png,patch2img} $(INSTALL)/bin
	@if test $(FORTRAN_UPPER) = 1; then \
	cp src/f90/include/*.MOD $(INSTALL)/include; \
	else \
	cp src/f90/include/*.mod $(INSTALL)/include; \
	fi
	@cp src/f90/scan_validate/scan_validate $(INSTALL)/bin
	@cp src/f90/l2gen/l2gen $(INSTALL)/bin
	@cp src/f90/l3gen/l3gen $(INSTALL)/bin
	@cp src/f90/ces_detect/ces_detect $(INSTALL)/bin

docs :
	@cd docs; $(MAKE)

libquiet :
	@cd src/f90/include; $(MAKE)

libquietscala:
	@cd src/f90/include; $(MAKE) libquietscala.a

tester :
	@cd src/f90/tester; $(MAKE)

postmap :
	@cd src/f90/postmap; $(MAKE)

map2cl :
	@cd src/f90/map2cl; $(MAKE)

scalapost :
	@cd src/f90/scalapost; $(MAKE)

map_editor :
	@cd src/f90/map_editor; $(MAKE)

tod2map :
	@cd src/f90/tod2map; $(MAKE)

scan_validate :
	@cd src/f90/scan_validate; $(MAKE)

l3gen :
	@cd src/f90/l3gen; $(MAKE)

scan_detect :
	@cd src/f90/scan_detect; $(MAKE)

l2gen :
	@cd src/f90/l2gen; $(MAKE)

utils: libutil
	@cd src/cpp/utils; $(MAKE)

libutil:
	@cd src/cpp/libutil; $(MAKE)

maptool:
	@cd src/f90/maptool; $(MAKE)

test :
	@cd src/f90/test; $(MAKE)

clean : clean_libquiet clean_postmap clean_map2cl clean_scalapost clean_map_editor clean_tod2map clean_libutil clean_utils clean_scan_validate clean_l3gen clean_scan_detect clean_l2gen clean_maptool clean_test

clean_postmap :
	@cd src/f90/postmap; $(MAKE) clean

clean_map2cl :
	@cd src/f90/map2cl; $(MAKE) clean

clean_scalapost :
	@cd src/f90/scalapost; $(MAKE) clean

clean_map_editor :
	@cd src/f90/map_editor; $(MAKE) clean

clean_tod2map :
	@cd src/f90/tod2map; $(MAKE) clean

clean_scan_validate :
	@cd src/f90/scan_validate; $(MAKE) clean 

clean_l3gen :
	@cd src/f90/l3gen; $(MAKE) clean 

clean_scan_detect :
	@cd src/f90/scan_detect; $(MAKE) clean 

clean_l2gen :
	@cd src/f90/l2gen; $(MAKE) clean 

clean_libquiet :
	@cd src/f90/include; $(MAKE) clean

clean_utils:
	@cd src/cpp/utils; $(MAKE) clean

clean_libutil :
	@cd src/cpp/libutil; $(MAKE) clean

clean_maptool :
	@cd src/f90/maptool; $(MAKE) clean

clean_tester :
	@cd src/f90/tester; $(MAKE) clean

clean_test :
	@cd src/f90/test; $(MAKE) clean

dist : clean
	@mkdir $(DIR)
	@mkdir -p $(DIR)/src/f90/include
	@mkdir -p $(DIR)/src/f90/postmap
	@mkdir -p $(DIR)/src/f90/map2cl
	@mkdir -p $(DIR)/src/f90/scalapost
	@mkdir -p $(DIR)/src/f90/map_editor
	@mkdir -p $(DIR)/src/f90/tod2map
	@mkdir -p $(DIR)/src/f90/scan_validate
	@mkdir -p $(DIR)/src/f90/l3gen
	@mkdir -p $(DIR)/src/f90/l2gen
	@mkdir -p $(DIR)/src/f90/scan_detect
	@mkdir -p $(DIR)/src/f90/utils
	@mkdir -p $(DIR)/src/cpp/utils
	@cp -r config Makefile $(DIR)
	@cp src/f90/include/*.f90 src/f90/include/Makefile $(DIR)/src/f90/include
	@cp src/f90/postmap/*.f90 src/f90/postmap/Makefile $(DIR)/src/f90/postmap
	@cp src/f90/map2cl/*.f90 src/f90/map2cl/Makefile $(DIR)/src/f90/map2cl
	@cp src/f90/scalapost/*.f90 src/f90/scalapost/Makefile $(DIR)/src/f90/scalapost
	@cp src/f90/map_editor/*.f90 src/f90/map_editor/Makefile $(DIR)/src/f90/map_editor
	@cp src/f90/tod2map/*.f90 src/f90/tod2map/Makefile $(DIR)/src/f90/tod2map
	@cp src/f90/scan_validate/*.f90 src/f90/scan_validate/Makefile $(DIR)/src/f90/scan_validate
	@cp src/f90/l3gen/*.f90 src/f90/l3gen/Makefile $(DIR)/src/f90/l3gen
	@cp src/f90/scan_detect/*.f90 src/f90/scan_detect/Makefile $(DIR)/src/f90/scan_detect
	@cp src/f90/l2gen/*.f90 src/f90/l2gen/Makefile $(DIR)/src/f90/l2gen
	@cp src/f90/utils/*.f90 src/f90/utils/Makefile $(DIR)/src/f90/utils
	@cp src/cpp/utils/{*.cpp,*.h,patch2img,Makefile} $(DIR)/src/cpp/utils
	@rm -rf $(DIR)/config/.svn
	@$(CTAR) $(DIR).tar.gz $(DIR)
	@rm -rf $(DIR)
