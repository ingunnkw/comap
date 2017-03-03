/*
   ephcom.h - header information for the ephemeris.com JPL Ephemeris software.
*/
#define EPHCOM_VERSION	"1.0"

#define EPHCOM_MAXLINE 128  /* Maximum # characters to allow in input line */
#define EPHCOM_MINJD -999999999.5
#define EPHCOM_MAXJD  999999999.5

/*
   Objects for pleph() ntarget and ncenter parameters.
*/
#define EPHCOM_MERCURY		1
#define EPHCOM_VENUS		2
#define EPHCOM_EARTH		3
#define EPHCOM_MARS		4
#define EPHCOM_JUPITER		5
#define EPHCOM_SATURN		6
#define EPHCOM_URANUS		7
#define EPHCOM_NEPTUNE		8
#define EPHCOM_PLUTO		9
#define EPHCOM_MOON		10 /* Moon, relative to Solar System center */
#define EPHCOM_SUN		11
#define EPHCOM_SSBARY		12 /* Solar System Barycenter */
#define EPHCOM_EMBARY		13 /* Earth-Moon Barycenter */
#define EPHCOM_NUTATION		14
#define EPHCOM_LIBRATION	15
#define EPHCOM_GEOMOON		16 /* Original Lunar Ephemeris coordinates */
#define EPHCOM_OBSERVER		17 /* User-defined observer position */
#define EPHCOM_NUMOBJECTS	17 /* Allocate memory for 17 solar sys objs */

/*
   This structure holds all the information contained in a JPLEPH header.
   When an ASCII or binary header is read, this structure is populated.
   Fill out this structure before writing an ASCII or binary header, and
   before performing any interpolations.
*/
struct ephcom_Header {
   int ksize;         /* block size, in first line of ASCII header */
   int ncoeff;        /* number of Chebyshev coefficients in data blocks */
   char ttl[3][86];   /* Hold up to 14*6=84 characters + "\n\0" */
   int ncon;          /* Number of assigned in cnam */
   char cnam[400][7]; /* Hold up to 400 6-character names ending with '\0' */
   int nval;          /* number of values for cval, to compare with ncon */
   double cval[400];  /* constant values, corresponding to cnam names */
   double au;         /* km/astronomical unit */
   double emrat;      /* Earth-Moon mass ratio */
   double clight;     /* Speed of light, km/sec */
   int numde;         /* ephemeris number */
   int numle;         /* lunar ephemeris number (can be same # as numde) */
   double ss[3];      /* start epoch, stop epoch, step size */
   int ipt[12][3];    /* index pointers into Chebyshev coefficients */
   int lpt[3];        /* libration pointer in a block */
   int maxcheby;      /* maximum Chebyshev coefficients for a body */
   };
/*
   This structure holds all interpolated positions of planets, Sun, and Moon
   at a given time.  All of the information available from interpolation
   at a given point in time is computed at once and preserved here.

   To populate this structure, have the ephemeris file open, and set:

              km, seconds, bary, et2[0], et2[1]

   Then call ephcom_get_coords() to get all coordinates, then call
   ephcom_pleph() for each desired (ntarget,ncenter) combination.
   See testeph.c for an example.  Note that unlike JPL's PLEPH()
   subroutine, you cannot call ephcom_pleph() without first initializing
   the pv[] array in ephcom_get_coords().

   There are some extra entries at the end, compared to JPL's PLEPH():

      pv[15][] - preserves the unmodified lunar ephemeris coordinates.
                 These coordinates are never offset by center, so that
                 their full precision is maintained.
      pv[16][] - start of possible user-defined coordinates.
                 These coordinates are offset by center.

   Note: None of the supplied routines read or modify the coordtype value.
*/
struct ephcom_Coords {
   int km;           /* 1 = positions in km; 0 = positions in AU            */
   int seconds;      /* 1 = timescale is seconds; 0 = timescale is days     */
   int bary;         /* 1 = Barycentric coordinates; 0 = adjust for center  */
   int center;       /* object to use as center (instead of SSBARY)         */
   int coordtype;    /* 0 = raw JPL ephemeris coordinates                   */
   double et2[2];    /* Julian Day of interpolation;
                        et2[0] = whole JD; et2[1] = fractional JD           */
   double pv[EPHCOM_NUMOBJECTS][6]; /* x, y, z Position & Velocity          */
                     /* pv[00..14][]: See Object numbers in #defines above */
                     /* pv[15][]: Geocentric Moon, from original lunar eph  */
                     /* pv[16][]: User-defined object                       */
   };

int ephcom_readbinary_header(FILE *, struct ephcom_Header *);
int ephcom_get_coords(FILE *, struct ephcom_Header *, struct ephcom_Coords *, double *);
int ephcom_pleph(struct ephcom_Coords *, int ntarg, int ncntr, double *);
