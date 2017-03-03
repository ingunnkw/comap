# We only depend on one library: fftw3

C := -I../libutil -I../ephcom -I. $(HEALPIX_INCLUDE_CPP) $(CFITSIO_INCLUDE) $(CFLAGS) $(FFTW_INCLUDE) $(HDF_INCLUDE)
L := -L../libutil -lutil -L../ephcom -lephcom $(HEALPIX_LINK_CPP) $(CFITSIO_LINK) $(FFTW_LINK) $(PNG_LINK) $(HDF_LINK_CPP)

targets = ppm2png powspec fft gp2tab mjd2unix tcat eph idate cut2accept ces_translate date2mjd l1_translate map2png #ps2fits ps2ascii patch_ident ces2runlist lxtxt ceslist_merge ceslistdiff osamize sfor hcat reverse ps2map map2ps ps2fits smap


all: $(targets)

%: %.cpp
	$(CXX) $C -o $@ $< $L

clean:
	rm -f $(targets) *.o
