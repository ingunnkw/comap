#
# Makefile for the ephemeris.com ephcom package to use JPL ephemeris files.
#
# This should run as is on Linux, and other Unix systems.  If there are
# any problems with a particular Unix system, email bugs@ephemeris.com.
#
# Set JPLDIR to the location of the JPL ephemeris files, copied from
# the CD from Willman-Bell or from JPL's FTP site:
#
#      ftp://ssd.jpl.nasa.gov/pub/eph/export
#
# You should be able to download and use the binary files from that site
# directly with this software, regardless of whether your host machine is
# Big-Endian or Little-Endian.
#
# This is a reduced version with only the ephcom library, for use in the
# Oslo QUIET pipeline.

CFLAGS = -O3

libephcom.a: gnulliver.o ephcom.o
	ar -rv $@ $^
	ranlib $@

clean:
	/bin/rm -rf *.o *.a
