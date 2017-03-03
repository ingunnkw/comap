# We only depend on one library: fftw3

C := -I../libutil -I../ephcom -I. $(HEALPIX_INCLUDE_CPP) $(CFITSIO_INCLUDE) $(CFLAGS) $(FFTW_INCLUDE) $(HDF_INCLUDE)
L := -L../libutil -lutil -L../ephcom -lephcom $(HEALPIX_LINK_CPP) $(CFITSIO_LINK) $(FFTW_LINK) $(PNG_LINK) $(HDF_LINK_CPP)

targets = map2png

all: $(targets)

%: %.cpp
	$(CXX) $C -o $@ $< $L

clean:
	rm -f $(targets) *.o
