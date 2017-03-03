/*
   gnulliver.c - Travels the middle of the road to swiftly gnullify
                 the Big-endian / Little-endian controversy.
                 Only works on machines with IEEE floating point.

   Swaps values between host order and network order without needing
   any include files and without calls to htonl(), etc.  This is a
   symmetric conversion, so the same function works in both directions.
   For example, to convert an integer between host byte ordering and
   network byte ordering (in either direction), use

         int i, j;
         j = gnulliver32(i);

   then use j in an equation (if you just read i from the net), or
   write j to the net (if i is an integer on your machine).

   There's no need to #include anything else -- just link and enjoy!

   The gnulliver routines work by setting a double to a certain value,
   then examining a unioned byte array.

   The functions are:

      char gnulliver() - determine (or recalculate) endianness of host:
         GNULLIVER_BIG, GNULLIVER_LITTLE, GNULLIVER_DEC
      unsigned short       gnulliver16(unsigned short input)
      unsigned             gnulliver32(unsigned input)
      unsigned long        gnulliver64(unsigned long input)
      unsigned long long   gnulliver128(unsigned long long input)
      float                gnulliver32f(float input)
      double               gnulliver64f(double input)
      long double          gnulliver128f(long double input)

   The code in the routines might appear largely redundant, and it is.
   Loops are also unrolled.  All this should minimize function calls and
   maximize speed.

   This library is distributed under the Lesser Gnu Public License.
   See http://www.gnu.org for a copy.

   Mail any bug reports, improvements, etc., to bugs@ephemeris.com
   Copyright 1994-2004 Paul Hardy.
*/
#define GNULLIVER_VERSION "1.0"

#define GNULLIVER_SWAP16  1   /* Swap adjacent bytes in 2-byte string */
#define GNULLIVER_SWAP32  2   /* Swap 2-byte pairs in 4-byte string   */
#define GNULLIVER_SWAP64  4   /* Swap 4-byte pairs in 8-byte string   */
#define GNULLIVER_SWAP128 8   /* Swap 8-byte pairs in 16-byte string  */
/*
   What we have to do to convert between network (Big-Endian) order and
   host machine order.
*/
#define GNULLIVER_BIG    0 /* No swapping if this machine is Big-endian */
#define GNULLIVER_LITTLE (GNULLIVER_SWAP16 | GNULLIVER_SWAP32 | GNULLIVER_SWAP64 | GNULLIVER_SWAP128)
/* This should work on a PDP-11 or VAX, but I have no way to test it. */
#define GNULLIVER_DEC    (GNULLIVER_SWAP32 | GNULLIVER_SWAP64 | GNULLIVER_SWAP128)
/*
   GNULLIVER_TEST is the magic number we use to determine byte swapping in
   a 64-bit word.  It assumes IEEE 754 floating point, and encodes to
   hexadecimal 0x40302010 00000000 on a Big-endian (network order) machine.
   The mantissa is large enough that it might work on a non-IEEE 754 compatible
   machine just by coincidence, but no guarantees.

   gnulliver currently assumes that if we swap 32-bit values in a 64-bit word,
   then we also swap 64-bit values in a 128-bit word (long double), which is
   true if the machine is Little-endian.  If this does not hold true on a
   particular platform, let me know and I'll modify the code.
*/
#define GNULLIVER_TEST (16.0 + (1.0/8.0 + 1.0/4096.0))

unsigned char gnulliver() {
   unsigned char result;
   union {
      double d;
      unsigned char ch[8];
      } bytes;
   int i;

   bytes.d = GNULLIVER_TEST;
   result = 0;

   if (bytes.ch[0] == 0) { /* swap 32-bit words in a 64-bit number */
      if (bytes.ch[4] == 0x10)               /* Assume 0x00000000 10203040 */
         result = GNULLIVER_LITTLE;
      else                                   /* Assume 0x00000000 20104030 */
         result = GNULLIVER_DEC;
      }
   else                                      /* Assume 0x40302010 00000000 */
      result = GNULLIVER_BIG;
// printf("Endian-ness: %x\n", result);
   return(result);
   }

unsigned short gnulliver16(unsigned short input) {
   static unsigned endian, init=1;
   unsigned char tmpch;
   union {
      unsigned short u;
      unsigned char ch[2];
      } bytes;

   if (init) {
      endian = gnulliver();
      init = 0;
      }
   if (endian != GNULLIVER_BIG)
      tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[1]; bytes.ch[1] = tmpch;

   return(bytes.u);
   }

unsigned gnulliver32(unsigned input) {
   static unsigned endian, init=1;
   char tmpch;
   union {
      unsigned u;
      unsigned char ch[4];
      } bytes;

   if (init) {
      endian = gnulliver();
      init = 0;
      }
   bytes.u = input;
   if (endian != GNULLIVER_BIG) {
      if (endian == GNULLIVER_LITTLE) {
         tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[3]; bytes.ch[3] = tmpch;
         tmpch = bytes.ch[1]; bytes.ch[1] = bytes.ch[2]; bytes.ch[2] = tmpch;
         }
      else { /* endian == GNULLIVER_DEC */
         tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[2]; bytes.ch[2] = tmpch;
         tmpch = bytes.ch[1]; bytes.ch[1] = bytes.ch[3]; bytes.ch[3] = tmpch;
         }
      }

   return(bytes.u);
   }

unsigned long gnulliver64(unsigned long input) {
   static unsigned endian, init=1;
   char tmpch;
   union {
      unsigned long u;
      char ch[8];
      } bytes;

   if (init) {
      endian = gnulliver();
      init = 0;
      }
   bytes.u = input;
   if (endian != GNULLIVER_BIG) {
      if (endian == GNULLIVER_LITTLE) {
         tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[7]; bytes.ch[7] = tmpch;
         tmpch = bytes.ch[1]; bytes.ch[1] = bytes.ch[6]; bytes.ch[6] = tmpch;
         tmpch = bytes.ch[2]; bytes.ch[2] = bytes.ch[5]; bytes.ch[5] = tmpch;
         tmpch = bytes.ch[3]; bytes.ch[3] = bytes.ch[4]; bytes.ch[4] = tmpch;
         }
      else { /* Assume GNULLIVER_DEC */
         tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[6]; bytes.ch[6] = tmpch;
         tmpch = bytes.ch[1]; bytes.ch[1] = bytes.ch[7]; bytes.ch[7] = tmpch;
         tmpch = bytes.ch[2]; bytes.ch[2] = bytes.ch[4]; bytes.ch[4] = tmpch;
         tmpch = bytes.ch[3]; bytes.ch[3] = bytes.ch[5]; bytes.ch[5] = tmpch;
         }
      }

   return(bytes.u);
   }


unsigned char *gnulliver64c(unsigned char *ch) {
   static unsigned endian, init=1;
   unsigned char tmpch;

   if (init) {
      endian = gnulliver();
      init = 0;
      }
// printf("\nG64c : %02x %02x %02x %02x %02x %02x %02x %02x ==>",
//        ch[0],ch[1],ch[2],ch[3],ch[4],ch[5],ch[6],ch[7]);
   if (endian != GNULLIVER_BIG) {
      if (endian == GNULLIVER_LITTLE) {
         tmpch = ch[0]; ch[0] = ch[7]; ch[7] = tmpch;
         tmpch = ch[1]; ch[1] = ch[6]; ch[6] = tmpch;
         tmpch = ch[2]; ch[2] = ch[5]; ch[5] = tmpch;
         tmpch = ch[3]; ch[3] = ch[4]; ch[4] = tmpch;
         }
      else { /* Assume GNULLIVER_DEC */
         tmpch = ch[0]; ch[0] = ch[6]; ch[6] = tmpch;
         tmpch = ch[1]; ch[1] = ch[7]; ch[7] = tmpch;
         tmpch = ch[2]; ch[2] = ch[4]; ch[4] = tmpch;
         tmpch = ch[3]; ch[3] = ch[5]; ch[5] = tmpch;
         }
      }
// printf("%02x %02x %02x %02x %02x %02x %02x %02x\n",
//        ch[0],ch[1],ch[2],ch[3],ch[4],ch[5],ch[6],ch[7]);

   return(ch);
   }

unsigned long long gnulliver128(unsigned long long input) {
   static unsigned endian, init=1;
   unsigned char tmpch;
   union {
      unsigned long long u;
      char ch[16];
      } bytes;

   if (init) {
      endian = gnulliver();
      init = 0;
      }
   bytes.u = input;
   if (endian != GNULLIVER_BIG) {
      if (endian == GNULLIVER_LITTLE) {
         tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[15]; bytes.ch[15] = tmpch;
         tmpch = bytes.ch[1]; bytes.ch[1] = bytes.ch[14]; bytes.ch[14] = tmpch;
         tmpch = bytes.ch[2]; bytes.ch[2] = bytes.ch[13]; bytes.ch[13] = tmpch;
         tmpch = bytes.ch[3]; bytes.ch[3] = bytes.ch[12]; bytes.ch[12] = tmpch;
         tmpch = bytes.ch[4]; bytes.ch[4] = bytes.ch[11]; bytes.ch[11] = tmpch;
         tmpch = bytes.ch[5]; bytes.ch[5] = bytes.ch[10]; bytes.ch[10] = tmpch;
         tmpch = bytes.ch[6]; bytes.ch[6] = bytes.ch[ 9]; bytes.ch[ 9] = tmpch;
         tmpch = bytes.ch[7]; bytes.ch[7] = bytes.ch[ 8]; bytes.ch[ 8] = tmpch;
         }
      else { /* Assume GNULLIVER_DEC */
         tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[14]; bytes.ch[14] = tmpch;
         tmpch = bytes.ch[1]; bytes.ch[1] = bytes.ch[15]; bytes.ch[15] = tmpch;
         tmpch = bytes.ch[2]; bytes.ch[2] = bytes.ch[12]; bytes.ch[12] = tmpch;
         tmpch = bytes.ch[3]; bytes.ch[3] = bytes.ch[13]; bytes.ch[13] = tmpch;
         tmpch = bytes.ch[4]; bytes.ch[4] = bytes.ch[10]; bytes.ch[10] = tmpch;
         tmpch = bytes.ch[5]; bytes.ch[5] = bytes.ch[11]; bytes.ch[11] = tmpch;
         tmpch = bytes.ch[6]; bytes.ch[6] = bytes.ch[ 8]; bytes.ch[ 8] = tmpch;
         tmpch = bytes.ch[7]; bytes.ch[7] = bytes.ch[ 9]; bytes.ch[ 9] = tmpch;
         }
      }

   return(bytes.u);
   }

float gnulliver32f(float input32) {
   union {
      float x32;
      unsigned char ch[4];
      } bytes;
   static unsigned endian, init=1;
   unsigned char tmpch;

   if (init) {
      endian = gnulliver();
      init = 0;
      }
   bytes.x32 = input32;
   if (endian != GNULLIVER_BIG) {
      if (endian == GNULLIVER_LITTLE) {
         tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[3]; bytes.ch[3] = tmpch;
         tmpch = bytes.ch[1]; bytes.ch[1] = bytes.ch[2]; bytes.ch[2] = tmpch;
         }
      else { /* Assume GNULLIVER_DEC */
         tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[2]; bytes.ch[2] = tmpch;
         tmpch = bytes.ch[1]; bytes.ch[1] = bytes.ch[3]; bytes.ch[3] = tmpch;
         }
      }

   return(bytes.x32);
   }

double gnulliver64f(double input64) {
   union {
      double x64;
      unsigned char ch[8];
      } bytes;
   static unsigned endian, init=1;
   unsigned char tmpch;
   int i;

   if (init) {
      endian = gnulliver();
      init = 0;
      }
   bytes.x64 = input64;
   if (endian != GNULLIVER_BIG) {
      if (endian == GNULLIVER_LITTLE) {
         tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[7]; bytes.ch[7] = tmpch;
         tmpch = bytes.ch[1]; bytes.ch[1] = bytes.ch[6]; bytes.ch[6] = tmpch;
         tmpch = bytes.ch[2]; bytes.ch[2] = bytes.ch[5]; bytes.ch[5] = tmpch;
         tmpch = bytes.ch[3]; bytes.ch[3] = bytes.ch[4]; bytes.ch[4] = tmpch;
         }
      else { /* Assume GNULLIVER_DEC */
         tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[6]; bytes.ch[6] = tmpch;
         tmpch = bytes.ch[1]; bytes.ch[1] = bytes.ch[7]; bytes.ch[7] = tmpch;
         tmpch = bytes.ch[2]; bytes.ch[2] = bytes.ch[4]; bytes.ch[4] = tmpch;
         tmpch = bytes.ch[3]; bytes.ch[3] = bytes.ch[5]; bytes.ch[5] = tmpch;
         }
      }

   return(bytes.x64);
   }

/*
   Note: This function's result will be unpredictable if the type
   long double is not a 128-bit word.
*/
long double gnulliver128f(long double input128) {
   union {
      long double x128;
      unsigned char ch[16];
      } bytes;
   static unsigned endian, init=1;
   unsigned char tmpch;

   if (init) {
      endian = gnulliver();
      init = 0;
      }
   bytes.x128 = input128;
   if (endian != GNULLIVER_BIG) {
      if (endian == GNULLIVER_LITTLE) {
         tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[15]; bytes.ch[15] = tmpch;
         tmpch = bytes.ch[1]; bytes.ch[1] = bytes.ch[14]; bytes.ch[14] = tmpch;
         tmpch = bytes.ch[2]; bytes.ch[2] = bytes.ch[13]; bytes.ch[13] = tmpch;
         tmpch = bytes.ch[3]; bytes.ch[3] = bytes.ch[12]; bytes.ch[12] = tmpch;
         tmpch = bytes.ch[4]; bytes.ch[4] = bytes.ch[11]; bytes.ch[11] = tmpch;
         tmpch = bytes.ch[5]; bytes.ch[5] = bytes.ch[10]; bytes.ch[10] = tmpch;
         tmpch = bytes.ch[6]; bytes.ch[6] = bytes.ch[ 9]; bytes.ch[ 9] = tmpch;
         tmpch = bytes.ch[7]; bytes.ch[7] = bytes.ch[ 8]; bytes.ch[ 8] = tmpch;
         }
      else { /* Assume GNULLIVER_DEC */
         tmpch = bytes.ch[0]; bytes.ch[0] = bytes.ch[14]; bytes.ch[14] = tmpch;
         tmpch = bytes.ch[1]; bytes.ch[1] = bytes.ch[15]; bytes.ch[15] = tmpch;
         tmpch = bytes.ch[2]; bytes.ch[2] = bytes.ch[12]; bytes.ch[12] = tmpch;
         tmpch = bytes.ch[3]; bytes.ch[3] = bytes.ch[13]; bytes.ch[13] = tmpch;
         tmpch = bytes.ch[4]; bytes.ch[4] = bytes.ch[10]; bytes.ch[10] = tmpch;
         tmpch = bytes.ch[5]; bytes.ch[5] = bytes.ch[11]; bytes.ch[11] = tmpch;
         tmpch = bytes.ch[6]; bytes.ch[6] = bytes.ch[ 8]; bytes.ch[ 8] = tmpch;
         tmpch = bytes.ch[7]; bytes.ch[7] = bytes.ch[ 9]; bytes.ch[ 9] = tmpch;
         }
      }

   return(bytes.x128);
   }
