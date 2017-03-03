========================================================
README file for the ephemeris.com JPL ephemeris software

                      Version 1.0
                      22 May 2004
========================================================


This package is a set of C routines and example programs
to access and use the data in the JPL ephemerides.  For
the latest version, check

     http://www.ephemeris.com/software/

The JPL DE405 ephemeris is available online at

     ftp://ssd.jpl.nasa.gov/pub/eph/export/

JPL ephemerides DE200, DE405, and DE406 are available on
CD from Willmann-Bell publishers, at

     http://www.willbell.com/

This software should run as is on any Unix/Linux or other
platform with Gnu's C compiler (gcc), as long as it uses IEEE
floating point representation for double precision numbers.
That is the floating point format of the binary ephemeris files.
If you're not sure about your machine, see the "make test" section
below.  If your machine uses a Pentium or Power PC processor,
it uses IEEE floating point.

This is free software; you can redistribute it and/or modify it under
the terms of version 2.1 of the GNU Lesser General Public License as
published by the Free Software Foundation.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MER-
CHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
General Public License for more details.  This is in the LICENSE.txt
file in this software distribution.

You should have received a copy of the GNU Lesser General Public
License along with this package; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
USA or see http://www.gnu.org/copyleft/lesser.html on the web.

If you do use this software in your own product, I'd appreciate your
mentioning that you use "the ephmeris.com software library" in your
documentation.

Please send mail to bugs@ephemeris.com to report any bugs.


=====================
SOFTWARE INSTALLATION
=====================
To compile all programs, just type

     make

There is no need to edit any files.  The programs should
run as is with any binary JPL ephemeris.

To both compile and install the software and Unix man pages:

1) Set the variable INSTALLDIR in file Makefile to
   the installation directory (by default, this is
   /usr/local/bin).
2) Set the variable MAN1DIR in file Makefile to the 
   directory where you want section 1 Unix man pages
   installed (by default ths is /usr/local/man/man1).
3) Type

     make install

4) To then remove object files (if you want), type

     make clean

To fully test software operation, you'll need the JPL
ephemeris CD from Willmann-Bell (and about 1.3 Gigabytes more
of free disk space for test output files).  Set the variable
JPLDIR to the directory containing the JPL ephemeris
files from the CD.  Then type

     make test

This will test operation with DE200, DE405, and DE406,
storing the results in a subdirectory named 'test'.
The Makefile creates this directory if it doesn't exist,
and removes files that were there if it does exist.
On a 1.7 GHz Pentium running Linux, this takes 11 minutes.

The tests convert all JPL binary ephemeris files into
ASCII files, back into binary files, and compare the results.
The tests also read the JPL TESTPO files to test for proper
interpolation.

If all you have is DE405 from the JPL FTP site, type

     make test405


====================
STAND ALONE PROGRAMS
====================
The following standalone programs are in this distribution.
They are compiled if you type 'make', and compiled and
installed along with man pages if you type 'make install'.
Consult the man pages or converted HTML pages for more
information on how to use these programs.

     eph2eph  - converts one or two binary ephemeris files
                into a third binary ephemeris file (this
                is similar to JPL's binshort.f and binmerge.f
                combined).
     eph2asc  - converts a binary ephemeris file into an ASCII
                header file and an ASCII data block file.
     asc2eph  - converts an ASCII header file and ASCII data
                block file to a binary ephemeris file.
     testeph  - interpolates positions in an ephemeris and
                compares results with JPL's TESTPO files from
                the JPL CD.
     headcmp  - prints two binary ephemeris headers formatted
                side by side to check for differences.
     ephcmp   - checks two binary ephemeris files for differences
                between constant values and data block values.
     ephcoeff - prints a binary ephemeris data block for a given
                Julian Day, printing a parsed, formatted output.
     vtransit - example program to adjust for Delta T and
                interpolate position of the 8 June 2004 Venus
                transit at mid-transit in rectangular coordinates.

The programs use three auxiliary files.

     gnulliver.c - routines to figure out the Big-endian/
                   Little-endian (and probably even VAX-endian)
                   byte order of a machine and swap bytes
                   to match network/Big-endian order (the
                   order of JPL's binary data).  These
                   routines are called within ephcom.c
                   library routines.  To write standalone
                   programs of your own, you won't need to
                   separately call any gnulliver routines.

     ephcom.h    - common definitions and structures.  Two
                   structures are defined here that are used
                   throughout the software:

                   ephcom_Header - data in an ephemeris header.
                   ephcom_Coords - interpolated positions of
                                   all bodies for a given time.

     ephcom.c    - the library of ephemeris routines (see the
                   ephcom.c file and example programs for more
                   on using these routines):

                   ephcom_readascii_header()   - read an ASCII header
                   ephcom_readascii_block()    - read an ASCII data block
                   ephcom_writeascii_header()  - write an ASCII header
                   ephcom_writeascii_block()   - write an ASCII data block
                   ephcom_readbinary_header()  - read a binary header
                   ephcom_readbinary_block()   - read a binary data block
                   ephcom_writebinary_header() - write a binary header
                   ephcom_writebinary_block()  - write a binary data block
                   ephcom_parse_block()        - print parsed binary data block
                   ephcom_nxtgrp()             - get next GROUP in ASCII header
                   ephcom_outdouble()          - write binary double precision
                   ephcom_outint()             - write binary integer
                   ephcom_indouble()           - read binary double precision
                   ephcom_inint()              - read binary integer
                   ephcom_doublstrc2f()        - convert C double format to FORTRAN
                   ephcom_pleph()              - compute object position & velocity
                   ephcom_get_coords()         - get all coordinates for a date
                   ephcom_cheby()              - perform Chebyshev interpolation
                   ephcom_jd2cal()             - convert Julian Day to calendar date
                   ephcom_cal2jd()             - convert calendar date to Julian Day


=============================
CHANGES FROM THE JPL SOFTWARE
=============================
If you've been using the JPL FORTRAN software, you will
notice a little similarity and a lot of difference between
this package and the JPL routines.  Some were necessary
by the nature of C, but most were design changes.

* COMMON blocks and many FORTRAN subroutine parameters
  have been converted to two C structures, ephcom_Header and
  ephcom_Coords, and they are defined in 'ephcom.h'.  There
  are no global variables.

* There is no FSIZER1(), FSIZER2(), or FSIZER3().  The ephcom.c
  routines read the ipt[][] and lpt[] index pointers in a header to
  determine the number of coefficients in an ephemeris, so you
  don't need to modify anything to read a new ephemeris file.

* There is no LIST array.  All positions and velocities for all
  objects on a given Julian Day are calculated at the same time.
  This takes less than 1 millisecond on a 1.5 GHz Pentium running
  Linux.

* There is no SPLIT().  Time in coords.et2[] is stored as
  coords.et2[0] = whole Julian Day, coords.et2[1] = fractional JD.

* There are no PVSUN or PNUT arrays.  Their coordinates are stored
  in pv[EPHCOM_SUN-1][] and pv[EPHCOM_NUTATION-1][], respectively.
  Libration is stored in pv[EPHCOM_LIBRATION-1][].

* CONST() has become ephcom_readascii_header() and
  ephcom_readbinary_header().

* There is no SELCON() subroutine; all constant names and values
  for a given header are stored in the header structure.

* STATE() and INTERP() have become ephcom_get_coords() and
  ephcom_cheby().  Whereas INTERP() is given an entire data block,
  ephcom_cheby() is only given coefficients for the current
  subinterval for the desired object.

* ephcom_cheby(), the Chebyshev interpolation routine, is given
  one double precision time normalized over the interval [-1,1],
  calculated from two double precision numbers stored in coords.et2[].
  It is not given two double precision numbers for time as in JPL's
  INTERP(), but the time only spans the desired subinterval, whereas
  with JPL's software the time is the entire Julian Day number.  For
  DE405 (with 32 day block spans), passing the entire date requires

       ceiling(log2(2500000 / 32)) = 17 bits more.

  For a 32-day subinterval using IEEE double precision floating point
  (1 sign bit + implied leading '1' bit + 52 mantissa bits), ephcom_cheby()
  time represents the limits of +/- 16 days (+/- 2^4 days) as:

       one sign bit
     + an implied leading '1' for the high order day bit
     + 3 remaining bits for the day, taken from the mantissa

  This leaves 52 - 3 = 49 mantissa bits for a fractional day.

       2^(-49) * 24 * 60 * 60 seconds = 0.153 nanoseconds.

  The worst-case time granularity of the interpolation should
  therefore be 0.153 nanoseconds.  If this does not suffice for
  some application, send a report to bugs@ephemeris.com.

* All coordinates (including those of the Moon) are converted to
  Solar System Barycentric coordinates, except of course nutations
  and librations.  The original geocentric Moon position and
  velocity is preserved in pv[EPHCOM_GEOMOON-1][] but never used
  by any of the example programs.  The 'testeph' program passes
  perfectly with the converted Solar System Barycentric Moon
  positions, but that only checks differences to 10^(-13).

* There are two unit conversion parameters for coordinates, not
  one (KM) as with JPL's software (note that the default is AU/day
  if both are zero):

     km=0: length is AU;          km=1: length is kilometers
     seconds=0: time is in days;  seconds=1: time is in seconds

  In JPL's software, setting KM to .TRUE. would convert length to
  kilometers and time to seconds.  The raw JPL ephemeris data itself
  is in kilometers/day (not kilometers/second).

* Before calling pleph() to compute rectangular coordinate
  distance and velocity between two objects, you must first call
  ephcom_get_coords(), which computes Solar System Barycentric
  positions and velocities at a given time.  This differs from
  JPL's PLEPH(), which itself directly calls STATE().

* FORTRAN double precision values are written in ASCII with an
  exponent indicator of 'D', not 'E'.  C uses 'E' for both single
  and double precision.  There is no way to reconcile these two.
  All ASCII ephemeris files write floating point numbers with a
  'D' exponent indicator, for compatibility with FORTRAN routines
  that might read the data.


===========================
JPL EPHEMERIS MODIFICATIONS
===========================
Some changes can be made to the original JPL header files.  All
changes were sanctioned by the creator of the JPL ephemerides,
Myles Standish.

* The last data block in a JPL binary ephemeris file is repeated
  once.  The ephemeris.com routines see that the time is the same
  as the next to last block, and ignore the final block.

* Most JPL ephemeris files give the ephemeris number as both the
  DE number and the Lunar Ephemeris (LE) number.  Some just use
  DE for both.  The ephemeris.com software writes the first title
  line as 'DEaaa/LEbbb', where aaa is taken from header.denum, and
  bbb is taken from header.lenum.

* The third header title line, which gives the ending Julian Day,
  ends exactly one block (32 or 64 days) after the actual ending
  day of the file.  The ephemeris.com routines that write binary
  and ASCII ephemeris headers regenerate the three title lines
  before writing them, by using the start and stop Julian Days
  in header.ss[0] and header.ss[1].  You'll notice this difference
  in the ASCII header files -- it is not a bug in the ephemeris.com
  routines.

* If there are no coefficients for an object in ipt[] or in lpt,
  then the first number in ipt or lpt for that object is set to the
  next available index location; it is not zero.  Some of the ephemeris
  files from JPL have zero in this location, and some have the next
  available index location.  Myles Standish indicated that the
  preferred format is to contain the next available location.

* Locations that aren't used at the end of the first and second block
  of a binary ephemeris file (the two header blocks) are filled with
  zeroes.  Some JPL binary ephemeris files contain non-zero data
  in these locations.

* The ephemeris.com eph2eph program can be used to remove the
  inconsistencies in an original JPL binary ephemeris file.  For
  example,

       eph2eph UNIX.405 UNIX.405 neweph.405

  will use JPL's UNIX.405 binary DE405 ephemeris from beginning to
  end as the single ephemeris file for input, and write a new ephemeris
  to neweph.405.


==========
WHAT NEXT?
==========
You've got a bunch of ultra-accurate rectangular coordinates -- now what?

More software exists, but hasn't been rigorously tested.  If you've
written software that uses the ephemeris.com package and would like
to contribute to it, send a report to bugs@ephemeris.com.  Also look
for new developments at

     http://www.ephemeris.com/software/

In the meantime, the U.S. Naval Observatory's "The Astronomical Almanac"
for the current year now gives instructions on converting rectangular
coordinates to topocentric coordinates, adjusted for light travel,
relativistic effects of the Sun, refraction, etc., along with tables of
positions for you to compare your results.  You can order the Almanac
online at http://www.amazon.com/.

Willmann-Bell also sells the book "Fundamental Ephemeris Computations"
that includes another (non-free) software package written in BASIC and C.
Their website is http://www.willbell.com/.
