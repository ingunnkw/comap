/*
   ephcom.c - software from http://www.ephemeris.com to
   access  and interpolate JPL Ephemeris data.
   Use this with the gnulliver.c routines for
   endian conversion.  Copyright 1994-2004 Paul Hardy.

   Note that this software is not a product of the Jet Propulsion
   Laboratory; it just uses and supports their ASCII and binary
   ephemeris files.  Please don't mail JPL concerning any bugs.
   Send bug reports or suggestions to bugs@ephemeris.com instead.

   The routines in this file and associated header files, accompanying
   programs, Makefile, and other files in the ephcom.tar/ephcom.zip package
   are released under the Lesser Gnu Public License (Lesser GPL).  This
   means that the code is free, with NO WARRANTY as to suitability for
   any platform expressed or implied.

   See http://www.gnu.org/copyleft/lesser.html for discussions
   on the implications of the Lesser GPL (LGPL).  See the end of
   this comment for the complete Lesser GPL.

   Because this is released under the Lesser GPL (LGPL), you can use
   this software in your own programs, even if you intend to sell those
   programs commercially.  I would appreciate your mentioning the
   "ephcom library from http://www.ephemeris.com" if you use these
   routines in your own software.

   Thanks, and enjoy!

   Paul Hardy, ephemeris.com, May 2004


   This file contains the following routines.  Running "make test" will
   check the proper operation of every one of these rouines.

   ephcom_readascii_header()   - read ASCII header file
   ephcom_readascii_block()    - read ASCII coefficient block
   ephcom_readbinary_header()  - read header from binary ephemeris
   ephcom_readbinary_block()   - read coefficient block from binary
   ephemeris
   ephcom_writeascii_header()  - write header in ASCII
   ephcom_writeascii_block()   - write coefficient block in ASCII
   ephcom_writebinary_header() - write header to binary ephemeris file
   ephcom_writebinary_block()  - write coefficient block to binary
   ephemeris file
   ephcom_parse_block()        - parse ("pretty print") a coefficient block
   ephcom_nxtgrp()             - read next "GROUP" from ASCII header
   ephcom_outdouble()          - write byte-swapped double to a file
   ephcom_outint()             - write byte-swapped int to a file
   ephcom_indouble()           - read byte-swapped double from a file
   ephcom_inint()              - read byte-swapped int from a file
   ephcom_doublstrc2f()        - change C ASCII double string to FORTRAN
   [there is no corresponding
   ephcom_doublestrf2c() routine;
   for FORTRAN to C format conversion,
   just change FORTRAN's double precision
   'D' exponent to 'E' in your software
   and everything else should parse fine]
   ephcom_pleph()              - calculate <x,y,z> and <xdot,ydot,zdot>
   for a given target and center, AFTER 
   calling ephcom_get_coords() (different
   sequence than with JPL's FORTRAN PLEPH)
   ephcom_get_coords()         - calculate <x,y,z> and <xdot,ydot,zdot>
   for all Solar System objects at a
   given time
   ephcom_cheby()              - interpolates Chebyshev coefficients
   for one sub-block of coefficients
   for one Solar System object at a
   given time
   ephcom_jd2cal()             - convert Julian Day to Julian or Gregorian
   Year, Month, Day, Hour, Minute, Second
ephcom_cal2jd()             - convert Julian or Gregorian calendar
Year, Month, Day, Hour, Minute, Second
to Julian Day

The indouble() and outdouble() routines rely upon the gnulliver64c()
	routine from gnulliver.c.

	The inint() and outint() routines rely upon the gnulliver32() routine
	from gnulliver.c.

	Below is the complete Lesser Gnu Public License Agreement, Version 2.1,
	under which all software in this package is distributed.


	GNU LESSER GENERAL PUBLIC LICENSE
	Version 2.1, February 1999


	Copyright (C) 1991, 1999 Free Software Foundation, Inc.
	59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
	Everyone is permitted to copy and distribute verbatim copies
	of this license document, but changing it is not allowed.

	[This is the first released version of the Lesser GPL.  It also counts
	as the successor of the GNU Library Public License, version 2, hence
	the version number 2.1.]

	Preamble
	The licenses for most software are designed to take away your freedom
	to share and change it. By contrast, the GNU General Public Licenses
	are intended to guarantee your freedom to share and change free
	software--to make sure the software is free for all its users.

	This license, the Lesser General Public License, applies to some
	specially designated software packages--typically libraries--of the
	Free Software Foundation and other authors who decide to use it. You
	can use it too, but we suggest you first think carefully about whether
	this license or the ordinary General Public License is the better
	strategy to use in any particular case, based on the explanations
	below.

	When we speak of free software, we are referring to freedom of use, not
	price. Our General Public Licenses are designed to make sure that you
	have the freedom to distribute copies of free software (and charge for
			this service if you wish); that you receive source code or can get it
	if you want it; that you can change the software and use pieces of it
	in new free programs; and that you are informed that you can do these
	things.

	To protect your rights, we need to make restrictions that forbid
	distributors to deny you these rights or to ask you to surrender these
	rights. These restrictions translate to certain responsibilities for
	you if you distribute copies of the library or if you modify it.

	For example, if you distribute copies of the library, whether gratis or
	for a fee, you must give the recipients all the rights that we gave
	you. You must make sure that they, too, receive or can get the source
	code. If you link other code with the library, you must provide
	complete object files to the recipients, so that they can relink them
	with the library after making changes to the library and recompiling
	it. And you must show them these terms so they know their rights.

	We protect your rights with a two-step method: (1) we copyright the
	library, and (2) we offer you this license, which gives you legal
	permission to copy, distribute and/or modify the library.

	To protect each distributor, we want to make it very clear that there
	is no warranty for the free library. Also, if the library is modified
	by someone else and passed on, the recipients should know that what
	they have is not the original version, so that the original author's
	reputation will not be affected by problems that might be introduced by
	others.

	Finally, software patents pose a constant threat to the existence of
	any free program. We wish to make sure that a company cannot
	effectively restrict the users of a free program by obtaining a
	restrictive license from a patent holder. Therefore, we insist that any
	patent license obtained for a version of the library must be consistent
	with the full freedom of use specified in this license.

	Most GNU software, including some libraries, is covered by the ordinary
	GNU General Public License. This license, the GNU Lesser General Public
	License, applies to certain designated libraries, and is quite
	different from the ordinary General Public License. We use this license
	for certain libraries in order to permit linking those libraries into
	non-free programs.

	When a program is linked with a library, whether statically or using a
	shared library, the combination of the two is legally speaking a
	combined work, a derivative of the original library. The ordinary
	General Public License therefore permits such linking only if the
	entire combination fits its criteria of freedom. The Lesser General
	Public License permits more lax criteria for linking other code with
	the library.

	We call this license the "Lesser" General Public License because it
	does Less to protect the user's freedom than the ordinary General
	Public License. It also provides other free software developers Less of
	an advantage over competing non-free programs. These disadvantages are
	the reason we use the ordinary General Public License for many
	libraries. However, the Lesser license provides advantages in certain
	special circumstances.

	For example, on rare occasions, there may be a special need to
	encourage the widest possible use of a certain library, so that it
	becomes a de-facto standard. To achieve this, non-free programs must be
	allowed to use the library. A more frequent case is that a free library
	does the same job as widely used non-free libraries. In this case,
	there is little to gain by limiting the free library to free software
	only, so we use the Lesser General Public License.

	In other cases, permission to use a particular library in non-free
	programs enables a greater number of people to use a large body of free
	software. For example, permission to use the GNU C Library in non-free
	programs enables many more people to use the whole GNU operating
	system, as well as its variant, the GNU/Linux operating system.

	Although the Lesser General Public License is Less protective of the
	users' freedom, it does ensure that the user of a program that is
	linked with the Library has the freedom and the wherewithal to run that
	program using a modified version of the Library.

	The precise terms and conditions for copying, distribution and
	modification follow. Pay close attention to the difference between a
	"work based on the library" and a "work that uses the library". The
	former contains code derived from the library, whereas the latter must
	be combined with the library in order to run.


	TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
	0. This License Agreement applies to any software library or other
	program which contains a notice placed by the copyright holder or other
	authorized party saying it may be distributed under the terms of this
	Lesser General Public License (also called "this License"). Each
	licensee is addressed as "you".

	A "library" means a collection of software functions and/or data
	prepared so as to be conveniently linked with application programs
	(which use some of those functions and data) to form executables.

	The "Library", below, refers to any such software library or work which
	has been distributed under these terms. A "work based on the Library"
	means either the Library or any derivative work under copyright law:
	that is to say, a work containing the Library or a portion of it,
	either verbatim or with modifications and/or translated
	straightforwardly into another language. (Hereinafter, translation is
			included without limitation in the term "modification".)

	"Source code" for a work means the preferred form of the work for
	making modifications to it. For a library, complete source code means
	all the source code for all modules it contains, plus any associated
	interface definition files, plus the scripts used to control
	compilation and installation of the library.

	Activities other than copying, distribution and modification are not
	covered by this License; they are outside its scope. The act of running
	a program using the Library is not restricted, and output from such a
	program is covered only if its contents constitute a work based on the
	Library (independent of the use of the Library in a tool for writing
			it). Whether that is true depends on what the Library does and what the
	program that uses the Library does.

	1. You may copy and distribute verbatim copies of the Library's
	complete source code as you receive it, in any medium, provided that
	you conspicuously and appropriately publish on each copy an appropriate
	copyright notice and disclaimer of warranty; keep intact all the
	notices that refer to this License and to the absence of any warranty;
	and distribute a copy of this License along with the Library.

	You may charge a fee for the physical act of transferring a copy, and
	you may at your option offer warranty protection in exchange for a fee.

	2. You may modify your copy or copies of the Library or any portion of
	it, thus forming a work based on the Library, and copy and distribute
	such modifications or work under the terms of Section 1 above, provided
	that you also meet all of these conditions:


	a) The modified work must itself be a software library.
	b) You must cause the files modified to carry prominent notices stating
	that you changed the files and the date of any change.
	c) You must cause the whole of the work to be licensed at no charge to
	all third parties under the terms of this License.
	d) If a facility in the modified Library refers to a function or a
	table of data to be supplied by an application program that uses the
	facility, other than as an argument passed when the facility is
	invoked, then you must make a good faith effort to ensure that, in the
	event an application does not supply such function or table, the
	facility still operates, and performs whatever part of its purpose
	remains meaningful.
	(For example, a function in a library to compute square roots has a
	 purpose that is entirely well-defined independent of the application.
	 Therefore, Subsection 2d requires that any application-supplied
	 function or table used by this function must be optional: if the
	 application does not supply it, the square root function must still
	 compute square roots.)

	These requirements apply to the modified work as a whole. If
	identifiable sections of that work are not derived from the Library,
	and can be reasonably considered independent and separate works in
	themselves, then this License, and its terms, do not apply to those
	sections when you distribute them as separate works. But when you
	distribute the same sections as part of a whole which is a work based
	on the Library, the distribution of the whole must be on the terms of
	this License, whose permissions for other licensees extend to the
	entire whole, and thus to each and every part regardless of who wrote
	it.

	Thus, it is not the intent of this section to claim rights or contest
	your rights to work written entirely by you; rather, the intent is to
	exercise the right to control the distribution of derivative or
	collective works based on the Library.

	In addition, mere aggregation of another work not based on the Library
	with the Library (or with a work based on the Library) on a volume of a
	storage or distribution medium does not bring the other work under the
	scope of this License.

	3. You may opt to apply the terms of the ordinary GNU General Public
	License instead of this License to a given copy of the Library. To do
	this, you must alter all the notices that refer to this License, so
	that they refer to the ordinary GNU General Public License, version 2,
	instead of to this License. (If a newer version than version 2 of the
			ordinary GNU General Public License has appeared, then you can specify
			that version instead if you wish.) Do not make any other change in
	these notices.

	Once this change is made in a given copy, it is irreversible for that
	copy, so the ordinary GNU General Public License applies to all
	subsequent copies and derivative works made from that copy.

	This option is useful when you wish to copy part of the code of the
	Library into a program that is not a library.

	4. You may copy and distribute the Library (or a portion or derivative
			of it, under Section 2) in object code or executable form under the
	terms of Sections 1 and 2 above provided that you accompany it with the
	complete corresponding machine-readable source code, which must be
	distributed under the terms of Sections 1 and 2 above on a medium
	customarily used for software interchange.

	If distribution of object code is made by offering access to copy from
	a designated place, then offering equivalent access to copy the source
	code from the same place satisfies the requirement to distribute the
	source code, even though third parties are not compelled to copy the
	source along with the object code.

	5. A program that contains no derivative of any portion of the Library,
	but is designed to work with the Library by being compiled or linked
	with it, is called a "work that uses the Library". Such a work, in
	isolation, is not a derivative work of the Library, and therefore falls
	outside the scope of this License.

	However, linking a "work that uses the Library" with the Library
	creates an executable that is a derivative of the Library (because it
			contains portions of the Library), rather than a "work that uses the
	library". The executable is therefore covered by this License. Section
	6 states terms for distribution of such executables.

	When a "work that uses the Library" uses material from a header file
	that is part of the Library, the object code for the work may be a
	derivative work of the Library even though the source code is not.
	Whether this is true is especially significant if the work can be
	linked without the Library, or if the work is itself a library. The
	threshold for this to be true is not precisely defined by law.

	If such an object file uses only numerical parameters, data structure
	layouts and accessors, and small macros and small inline functions (ten
			lines or less in length), then the use of the object file is
	unrestricted, regardless of whether it is legally a derivative work.
	(Executables containing this object code plus portions of the Library
	 will still fall under Section 6.)

	Otherwise, if the work is a derivative of the Library, you may
	distribute the object code for the work under the terms of Section 6.
	Any executables containing that work also fall under Section 6, whether
	or not they are linked directly with the Library itself.

	6. As an exception to the Sections above, you may also combine or link
	a "work that uses the Library" with the Library to produce a work
	containing portions of the Library, and distribute that work under
	terms of your choice, provided that the terms permit modification of
	the work for the customer's own use and reverse engineering for
	debugging such modifications.

	You must give prominent notice with each copy of the work that the
	Library is used in it and that the Library and its use are covered by
	this License. You must supply a copy of this License. If the work
	during execution displays copyright notices, you must include the
	copyright notice for the Library among them, as well as a reference
	directing the user to the copy of this License. Also, you must do one
	of these things:


	a) Accompany the work with the complete corresponding machine-readable
	source code for the Library including whatever changes were used in the
	work (which must be distributed under Sections 1 and 2 above); and, if
	the work is an executable linked with the Library, with the complete
	machine-readable "work that uses the Library", as object code and/or
	source code, so that the user can modify the Library and then relink to
	produce a modified executable containing the modified Library. (It is
			understood that the user who changes the contents of definitions files
			in the Library will not necessarily be able to recompile the
			application to use the modified definitions.)
	b) Use a suitable shared library mechanism for linking with the
	Library. A suitable mechanism is one that (1) uses at run time a copy
	of the library already present on the user's computer system, rather
	than copying library functions into the executable, and (2) will
	operate properly with a modified version of the library, if the user
	installs one, as long as the modified version is interface-compatible
	with the version that the work was made with.
	c) Accompany the work with a written offer, valid for at least three
	years, to give the same user the materials specified in Subsection 6a,
	above, for a charge no more than the cost of performing this
	distribution.
	d) If distribution of the work is made by offering access to copy from
	a designated place, offer equivalent access to copy the above specified
	materials from the same place.
	e) Verify that the user has already received a copy of these materials
	or that you have already sent this user a copy.
	For an executable, the required form of the "work that uses the
	Library" must include any data and utility programs needed for
	reproducing the executable from it. However, as a special exception,
	the materials to be distributed need not include anything that is
	normally distributed (in either source or binary form) with the major
	components (compiler, kernel, and so on) of the operating system on
	which the executable runs, unless that component itself accompanies the
	executable.

	It may happen that this requirement contradicts the license
	restrictions of other proprietary libraries that do not normally
	accompany the operating system. Such a contradiction means you cannot
	use both them and the Library together in an executable that you
	distribute.

	7. You may place library facilities that are a work based on the
	Library side-by-side in a single library together with other library
	facilities not covered by this License, and distribute such a combined
	library, provided that the separate distribution of the work based on
	the Library and of the other library facilities is otherwise permitted,
	and provided that you do these two things:


	a) Accompany the combined library with a copy of the same work based on
	the Library, uncombined with any other library facilities. This must be
	distributed under the terms of the Sections above.
	b) Give prominent notice with the combined library of the fact that
	part of it is a work based on the Library, and explaining where to find
	the accompanying uncombined form of the same work.
	8. You may not copy, modify, sublicense, link with, or distribute the
	Library except as expressly provided under this License. Any attempt
	otherwise to copy, modify, sublicense, link with, or distribute the
	Library is void, and will automatically terminate your rights under
	this License. However, parties who have received copies, or rights,
	from you under this License will not have their licenses terminated so
	long as such parties remain in full compliance.

	9. You are not required to accept this License, since you have not
	signed it. However, nothing else grants you permission to modify or
	distribute the Library or its derivative works. These actions are
	prohibited by law if you do not accept this License. Therefore, by
	modifying or distributing the Library (or any work based on the
			Library), you indicate your acceptance of this License to do so, and
	all its terms and conditions for copying, distributing or modifying the
	Library or works based on it.

	10. Each time you redistribute the Library (or any work based on the
			Library), the recipient automatically receives a license from the
	original licensor to copy, distribute, link with or modify the Library
	subject to these terms and conditions. You may not impose any further
	restrictions on the recipients' exercise of the rights granted herein.
	You are not responsible for enforcing compliance by third parties with
	this License.

	11. If, as a consequence of a court judgment or allegation of patent
	infringement or for any other reason (not limited to patent issues),
	conditions are imposed on you (whether by court order, agreement or
			otherwise) that contradict the conditions of this License, they do not
	excuse you from the conditions of this License. If you cannot
	distribute so as to satisfy simultaneously your obligations under this
	License and any other pertinent obligations, then as a consequence you
	may not distribute the Library at all. For example, if a patent license
	would not permit royalty-free redistribution of the Library by all
	those who receive copies directly or indirectly through you, then the
	only way you could satisfy both it and this License would be to refrain
	entirely from distribution of the Library.

	If any portion of this section is held invalid or unenforceable under
	any particular circumstance, the balance of the section is intended to
	apply, and the section as a whole is intended to apply in other
	circumstances.

	It is not the purpose of this section to induce you to infringe any
	patents or other property right claims or to contest validity of any
	such claims; this section has the sole purpose of protecting the
	integrity of the free software distribution system which is implemented
	by public license practices. Many people have made generous
	contributions to the wide range of software distributed through that
	system in reliance on consistent application of that system; it is up
	to the author/donor to decide if he or she is willing to distribute
	software through any other system and a licensee cannot impose that
	choice.

	This section is intended to make thoroughly clear what is believed to
	be a consequence of the rest of this License.

	12. If the distribution and/or use of the Library is restricted in
	certain countries either by patents or by copyrighted interfaces, the
	original copyright holder who places the Library under this License may
	add an explicit geographical distribution limitation excluding those
	countries, so that distribution is permitted only in or among countries
	not thus excluded. In such case, this License incorporates the
	limitation as if written in the body of this License.

	13. The Free Software Foundation may publish revised and/or new
	versions of the Lesser General Public License from time to time. Such
	new versions will be similar in spirit to the present version, but may
	differ in detail to address new problems or concerns.

	Each version is given a distinguishing version number. If the Library
	specifies a version number of this License which applies to it and "any
	later version", you have the option of following the terms and
	conditions either of that version or of any later version published by
	the Free Software Foundation. If the Library does not specify a license
	version number, you may choose any version ever published by the Free
	Software Foundation.

	14. If you wish to incorporate parts of the Library into other free
	programs whose distribution conditions are incompatible with these,
	write to the author to ask for permission. For software which is
	copyrighted by the Free Software Foundation, write to the Free Software
	Foundation; we sometimes make exceptions for this. Our decision will be
	guided by the two goals of preserving the free status of all
	derivatives of our free software and of promoting the sharing and reuse
	of software generally.

	NO WARRANTY

	15. BECAUSE THE LIBRARY IS LICENSED FREE OF CHARGE, THERE IS NO
	WARRANTY FOR THE LIBRARY, TO THE EXTENT PERMITTED BY APPLICABLE LAW.
	EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR
	OTHER PARTIES PROVIDE THE LIBRARY "AS IS" WITHOUT WARRANTY OF ANY KIND,
	EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
	ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE LIBRARY IS WITH
	YOU. SHOULD THE LIBRARY PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
	NECESSARY SERVICING, REPAIR OR CORRECTION.

	16. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN
	WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY
	AND/OR REDISTRIBUTE THE LIBRARY AS PERMITTED ABOVE, BE LIABLE TO YOU
	FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR
	CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE
	LIBRARY (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
			RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
			FAILURE OF THE LIBRARY TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
	SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
	DAMAGES.


	END OF TERMS AND CONDITIONS

	*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "ephcom.h"

	int ephcom_outdouble(FILE *, double);

	/*
	   Read a JPL Ephemeris ASCII header from the file pointed to by infp
	   and store values in header structure.  Write any errors to stderr.
	 */
	ephcom_readascii_header(FILE *infp, struct ephcom_Header *header) {

		char group[13]; /* To store the "GROUP" header line information */
		double val1, val2, val3; /* To read text line with 3 double precision words */
		int i, j, k, n;
		int iword; /* word number we're reading in a line */
		int blockout; /* number of bytes we've written to current block/rec in file */
		int blockbytes; /* number of bytes in a block, equals 8 * ncoeff */

		char readbuf[EPHCOM_MAXLINE + 1];

		char *ephcom_nxtgrp(char *, char *, FILE *);
		char *fgets(char *, int, FILE *);

		char outhcars[EPHCOM_MAXLINE + 1];
		size_t fwrite(const void *ptr, size_t size, size_t  nmemb, FILE *stream);

		/*
		   First header line: KSIZE= # NCOEFF= #
		 */
		if (infp != stdin) rewind(infp);
		fgets(readbuf, EPHCOM_MAXLINE, infp);
		sscanf(readbuf, "%*6s%6d%*11s%6d", &header->ksize, &header->ncoeff);
		blockbytes = 8 * header->ncoeff;  /* The size of a double, times # of doubles/block */
		if (header->ksize != 2*header->ncoeff) {
			fprintf(stderr, "Badly formed header; KSIZE <> 2*NCOEFF\n\n");
			exit(1);
		}
		/*
		   GROUP 1010: Title of ephemeris (DE/LE number, start JD, end JD)
		 */
		ephcom_nxtgrp(group, "GROUP   1010", infp);
		fgets(header->ttl[0], EPHCOM_MAXLINE, infp);  /* JPL Ephemeris title line */
		if (strncmp(header->ttl[0], "JPL ", 4) != 0) {
			fprintf(stderr,"\nERROR: file is not a JPL ASCII header file.\n\n");
			exit(1);
		}
		fgets(header->ttl[1], EPHCOM_MAXLINE, infp);  /* Start epoch */
		fgets(header->ttl[2], EPHCOM_MAXLINE, infp);  /* Finish epoch */
		/*
		   Convert any newlines or tabs to single spaces.
		 */
		for (i=0; i<3; i++) {
			for (j=0; j<84; j++)
				if (isspace(header->ttl[i][j]))  header->ttl[i][j] = ' ';
			header->ttl[i][84] = '\0';
		}
		/*
		   GROUP 1030: Start and End JD, timestep (in JD) per block.
		 */
		ephcom_nxtgrp(group, "GROUP   1030", infp);
		fgets(readbuf, EPHCOM_MAXLINE, infp);
		sscanf(readbuf, " %lE %lE %lE", &header->ss[0], &header->ss[1], &header->ss[2]);
		/*
		   GROUP 1040: Constant names.
		 */
		ephcom_nxtgrp(group, "GROUP   1040", infp);
		fgets(readbuf, EPHCOM_MAXLINE, infp);
		header->ncon = atoi(readbuf);
		/*
		   Now read the constant names, 10 per line, each 6 characters long
		   preceded by 2 blanks.  Pad names with blanks to make 6 characters.
		 */
		for (i=0; i<header->ncon; ) {
			fgets(readbuf, EPHCOM_MAXLINE, infp);
			if ((j = strlen(readbuf)) < 81) {  /* Pad end with blanks for copying */
				while (j < 81) readbuf[j] = ' ';
			}
			for (iword=0; iword<10 && i<header->ncon; iword++, i++) {
				strncpy(header->cnam[i], &readbuf[2 + iword*8], 6);
				header->cnam[i][6] = '\0';
			}
		}
		/*
		   GROUP 1041: Constant values.
		 */
		ephcom_nxtgrp(group, "GROUP   1041", infp);
		fgets(readbuf, EPHCOM_MAXLINE, infp);
		header->nval = atoi(readbuf);
		if (header->nval != header->ncon) {
			fprintf(stderr,"Error: number of constants and values not equal.\n\n");
			exit(1);
		}
		/*
		   Now read constant values, 3 per line, 26 characters each.
		 */
		for (i = 0; i < header->ncon; i += 3) {
			fgets(readbuf, EPHCOM_MAXLINE, infp);
			for (j=0; j<strlen(readbuf); j++)
				if (tolower(readbuf[j]) == 'd') readbuf[j] = 'E'; /* exponent is 'E' */
			sscanf(readbuf, "%lE %lE %lE",
					&header->cval[i], &header->cval[i+1], &header->cval[i+2]);
		}
		/*
		   GROUP 1050: Constant values.
		 */
		ephcom_nxtgrp(group, "GROUP   1050", infp);
		for (i =0; i < 3; i++) {
			fgets(readbuf, EPHCOM_MAXLINE, infp); /* Read line of 13 6-digit integers */
			for (j = 0; j < 12; j++) {
				header->ipt[j][i] = atoi(&readbuf[6*j]);
			}
			header->lpt[i] = atoi(&readbuf[6*12]);
		}
		/*
		   If there are no coefficients for an ipt[i][] object (i.e., ipt[i][1]==0),
		   then ipt[i][0] should contain the value of the next available coefficient
		   number rather than 0, as per communication of Myles Standish to Paul Hardy
		   on preferred format of ephemeris headers.

		   If there are no libration coefficients (i.e., lpt[1]==0), then lpt[0]
		   should contain the value of the next available coefficient number rather
		   than 0 as well, as per the same communication from Myles Standish.
		 */
		/* First set j to maximum index into ipt[] that has coefficients */
		j = 0;
		for (i=1; i<12; i++)
			if (header->ipt[i][1] > 0 && header->ipt[i][0] > j)
				j = i;
		/* Now set j to next available index count. */
		if (header->lpt[1] > 0 && header->lpt[0] > j)
			j = header->lpt[1] + header->lpt[1] * header->lpt[2] * 3;
		else
			j = header->ipt[j][0] +
				header->ipt[j][1] * header->ipt[j][2] * (j==11 ? 2 : 3);
		for (i=1; i<12; i++)
			if (header->ipt[i][0] == 0) header->ipt[i][0] = j;
		if (header->lpt[0] == 0) header->lpt[0] = j;
		/*
		   Set the maximum number of Chebyshev coefficients possible for this file,
		   to initialize position and velocity Chebyshev coefficient arrays during
		   Chebyshev interpolation.
		 */
		header->maxcheby = 0;
		for (i=0; i<12; i++)
			if (header->ipt[i][1] > header->maxcheby)
				header->maxcheby = header->ipt[i][1];
		if (header->lpt[1] > header->maxcheby)
			header->maxcheby = header->lpt[1];

		header->au = 0.0;
		header->emrat = 0.0;
		header->numde = 0;
		for (i = 0; i < header->ncon; i++) {
			if (strncmp(header->cnam[i], "AU    ", 6) == 0)
				header->au = header->cval[i];
			else if (strncmp(header->cnam[i], "EMRAT ", 6) == 0)
				header->emrat = header->cval[i];
			else if (strncmp(header->cnam[i], "DENUM ", 6) == 0)
				header->numde = header->cval[i];
			else if (strncmp(header->cnam[i], "CLIGHT", 6) == 0)
				header->clight = header->cval[i];
			else if (strncmp(header->cnam[i], "LENUM ", 6) == 0)
				header->numle = header->cval[i];
		}
		if (header->numle == 0) header->numle = header->numde;
		/*
		   GROUP 1070: Constant values.
		 */
		ephcom_nxtgrp(group, "GROUP   1070", infp);
		/*
		   Now we're pointing to the first block of coefficient data, after header.
		   Return at the point where we can start reading coefficients.
		 */
		return(0);
	}

/*
   Read a block of data coefficients from a JPL ASCII ephemeris file.
   Returns number of coefficients read, 0 at EOF.
 */
ephcom_readascii_block(
		FILE *infp,
		struct ephcom_Header *header,
		double *datablock) {

	int i,j;
	int datalines; /* lines of data we've read */
	int datapoints; /* points of data we've read/converted/written */
	char readbuf[EPHCOM_MAXLINE + 1];
	double val1, val2, val3; /* To read text line with 3 double precision words */

	/*
	   First line in an ASCII block will be the block number, followed by
	   the number of coefficients.
	 */
	datalines = 0;  /* Not reported, but leave in for debugging */
	datapoints = 0;
	if (fgets(readbuf, EPHCOM_MAXLINE, infp) && !feof(infp)) {
		sscanf(readbuf, "%d %d", &i, &j);
		if (j != header->ncoeff) {
			fprintf(stderr,
					"\nERROR: ASCII data file's %d coefficients/block\n", j);
			fprintf(stderr,
					"       doesn't match ASCII header's %d coefficients/block.\n\n",
					header->ncoeff);
			exit(1);
		}
		datalines++;
		for (i=0; i < header->ncoeff && !feof(infp); i += 3) {
			fgets(readbuf, EPHCOM_MAXLINE, infp);
			for (j=0; j<strlen(readbuf); j++)
				if (tolower(readbuf[j]) == 'd')  readbuf[j] = 'e';
			datalines++;
			/*
			   This is horrible, but use "%le" here and "%lE in the other
			   ASCII data routine (ephcom_readascii_header) so gcc won't try
			   to store the formats in the same location and write to them.
			   (Problem with gcc not acting like K&R without -traditional flag
			   and without -fwritable-strings flag.)
			 */
			sscanf(readbuf, " %le %le %le", &val1, &val2, &val3);
			datablock[i] = val1;
			datapoints++;
			if ((i+1) < header->ncoeff) {
				datablock[i+1] = val2;
				datapoints++;
				if ((i+2) < header->ncoeff) {
					datablock[i+2] = val3;
					datapoints++;
				}
			}
		}
	}
	return(datapoints);
}

/*
   Read a JPL Ephemeris header in binary format.  Store values in
   an ephcom_Header struct.
 */
ephcom_readbinary_header(FILE *infp, struct ephcom_Header *header) {

	int i, j, k;
	char *fgets(char *, int, FILE *);
	int fseek(FILE *, long, int);
	int fgetc(FILE *);
	double ephcom_indouble(FILE *);
	int ephcom_inint(FILE *);

	if (infp != stdin) rewind(infp);
	/*
	   Read title lines.
	 */
	for (i=0; i<3; i++) {
		for (j=0; j<84; j++) {
			header->ttl[i][j] = fgetc(infp);
		}
		if (i == 0 && strncmp(header->ttl[0], "JPL ", 4) != 0) {
			fprintf(stderr,"\nERROR: file is not a JPL ephemeris file.\n\n");
			if (strncmp(header->ttl[0], "KSIZE", 5) == 0)
				fprintf(stderr,"File is an ASCII JPL ephemeris header instead.\n\n");
			exit(1);
		}
		header->ttl[i][j] = '\0';
	}
	/*
	   Read constant names.
	 */
	for (i=0; i<400; i++) {
		for (j=0; j<6; j++) {
			header->cnam[i][j] = fgetc(infp);
		}
		header->cnam[i][j] = '\0';
	}
	/*
	   Read ephemeris start epoch, stop epoch, and step size (in Julian Days).
	 */
	for (i=0; i<3; i++) {
		header->ss[i] = ephcom_indouble(infp);
	}
	/*
	   Read NCON, AU, EMRAT.
	 */
	header->ncon  = ephcom_inint(infp);
	header->au    = ephcom_indouble(infp);
	header->emrat = ephcom_indouble(infp);
	header->nval  = header->ncon;
	/*
	   Read indexes for coefficients in data block.  Written in transposed
	   order (Fortran and C matrices are transposed).
	 */
	for (i=0; i<12; i++) {
		for (j=0; j<3; j++) {
			header->ipt[i][j] = ephcom_inint(infp);
		}
	}
	header->numde = ephcom_inint(infp);  /* Get ephemeris number */
	for (i=0; i<3; i++) header->lpt[i] = ephcom_inint(infp);
	/*
	   If there are no coefficients for an ipt[i][] object (i.e., ipt[i][1]==0),
	   then ipt[i][0] should contain the value of the next available coefficient
	   number rather than 0, as per communication of Myles Standish to Paul Hardy
	   on preferred format of ephemeris headers.

	   If there are no libration coefficients (i.e., lpt[1]==0), then lpt[0]
	   should contain the value of the next available coefficient number rather
	   than 0 as well, as per the same communication from Myles Standish.
	 */
	/* First set j to maximum index into ipt[] that has coefficients */
	j = 0;
	for (i=1; i<12; i++)
		if (header->ipt[i][1] > 0 && header->ipt[i][0] > j)
			j = i;
	/* Now set j to next available index count. */
	if (header->lpt[1] > 0 && header->lpt[0] > j)
		j = header->lpt[1] + header->lpt[1] * header->lpt[2] * 3;
	else
		j = header->ipt[j][0] +
			header->ipt[j][1] * header->ipt[j][2] * (j==11 ? 2 : 3);
	for (i=1; i<12; i++)
		if (header->ipt[i][0] == 0) header->ipt[i][0] = j;
	if (header->lpt[0] == 0) header->lpt[0] = j;
	/*
	   Set the maximum number of Chebyshev coefficients possible for this file,
	   to initialize position and velocity Chebyshev coefficient arrays during
	   Chebyshev interpolation.
	 */
	header->maxcheby = 0;
	for (i=0; i<12; i++)
		if (header->ipt[i][1] > header->maxcheby)
			header->maxcheby = header->ipt[i][1];
	if (header->lpt[1] > header->maxcheby)
		header->maxcheby = header->lpt[1];

	// /*
	//    From JPL ephemeris number, set NCOEFF and calculate KSIZE = 2*NCOEFF.
	// */
	// switch (header->numde) {
	//    case 102:
	//       header->ncoeff = 773;
	//       break;
	//    case 200:
	//       header->ncoeff = 826;
	//       break;
	//    case 202:
	//       header->ncoeff = 826;
	//       break;
	//    case 403:
	//       header->ncoeff = 1018;
	//       break;
	//    case 405:
	//       header->ncoeff = 1018;
	//       break;
	//    case 406:
	//       header->ncoeff = 728;
	//       break;
	//    case 410:
	//       header->ncoeff = 1018;
	//       break;
	//    default:
	//       header->ncoeff = 1018;
	//       break;
	//    }
	/*
	   Calculate number of coefficients, starting with
	   highest index into a data block (stored in j).
	 */
	j = 0;
	for (i=1; i<12; i++)
		if (header->ipt[i][1] > 0 && header->ipt[i][0] > header->ipt[j][0]) j = i;

	/*
	   Now see if the starting point we found is lower than where
	   lpt[] starts.  If not, use lpt[] for largest value.
	 */
	if (header->lpt[1] > 0 && header->lpt[0] > header->ipt[j][0]) {
		header->ncoeff = header->lpt[0] - 1 +     /* starting point */
			(header->lpt[1] *        /* coefficients per coordinate */
			 header->lpt[2]) *       /* subblocks per block */
			3;                      /* coordinates */
	}
	else {
		header->ncoeff = header->ipt[j][0] - 1 +  /* starting point */
			(header->ipt[j][1] *     /* coefficients per coordinate */
			 header->ipt[j][2]) *    /* subblocks per block */
			(j == 11 ? 2 : 3);       /* coordinates */
	}
	// printf("Number of coefficients: %d\n", header->ncoeff);

	header->ksize = header->ncoeff + header->ncoeff; /* KSIZE = 2*NCOEFF always */
	/*
	   Skip to second block in file.
	 */
	fseek(infp, header->ncoeff * 8, SEEK_SET);
	/*
	   Read ephemeris constants.
	 */
	for (i=0; i<header->ncon; i++) {
		header->cval[i] = ephcom_indouble(infp);
		if (strncmp(header->cnam[i], "NCOEFF", 6) == 0) {
			header->ncoeff = header->cval[i];
			header->ksize  = 2 * header->cval[i];
		}
		else if (strncmp(header->cnam[i], "LENUM ", 6) == 0)
			header->numle = header->cval[i];
	}
	if (header->numle == 0) header->numle = header->numde;

	return(0);
}

/*
   Read a JPL Ephemeris data block in binary format.

   This is the only routine in this library that accesses a file
   as a direct access file, with a specified block number.  The
   block number ranges from 0 on up (starting at first data block,
   after the 2 header blocks).  Returns the number of coefficients
   read, or 0 at EOF.
 */
int ephcom_readbinary_block(
		FILE *infp,                    /* File pointer for direct access file */
		struct ephcom_Header *header,  /* header struct, already filled in    */
		int blocknum,                  /* Data block number, starting with 0  */
		double *datablock              /* returned coefficient data block     */
		) {

	int i;
	long filebyte;
	double ephcom_indouble(FILE *);
	int fseek(FILE *, long, int);

	filebyte = (blocknum + 2) * header->ncoeff * 8; /* 8 bytes per coefficient */
	fseek(infp, filebyte, SEEK_SET);
	// printf("Blocknum: %d, Byte: %d\n", blocknum, filebyte);
	for (i=0; !feof(infp) && i<header->ncoeff; i++) {
		datablock[i] = ephcom_indouble(infp);
	}
	if (i < header->ncoeff && feof(infp)) i = -1; /* 0 --> EOF */
	return(i); /* Number of coefficients successfuly read (all or nohing). */
}


/*
   Write header information in ASCII format.
 */
ephcom_writeascii_header(FILE *outfp, struct ephcom_Header *header) {

	char group[13];
	double val1, val2, val3; /* To read text line with 3 double precision words */
	int i, j, k, n;
	int iword; /* word number we're reading in a line */
	int blockout; /* number of bytes we've written to current block/rec in file */
	int blockbytes; /* number of bytes in a block, equals 8 * ncoeff */
	static char spaces[84]="                                                                                \r\n";
	int idate[6];
	char *month[12] = {"JAN","FEB","MAR","APR","MAY","JUN",
		"JUL","AUG","SEP","OCT","NOV","DEC"};

	char writebuf[EPHCOM_MAXLINE + 1];
	char outhcars[EPHCOM_MAXLINE + 1];
	size_t fwrite(const void *ptr, size_t size, size_t  nmemb, FILE *stream);

	/*
	   First header line: KSIZE= # NCOEFF= #
	 */
	blockbytes = 8 * header->ncoeff;  /* sizeof(double) * # of doubles/block */
	sprintf(writebuf, "KSIZE=%6d    NCOEFF=%6d", header->ksize, header->ncoeff);
	k = strlen(writebuf);
	strcpy(&writebuf[k], &spaces[k]);
	fprintf(outfp, writebuf);
	if (header->ksize != 2*header->ncoeff) {
		fprintf(stderr, "Badly formed header; KSIZE <> 2*NCOEFF\r\n\r\n");
		exit(1);
	}
	/*
	   GROUP 1010: Title of ephemeris (DE/LE number, start JD, end JD)
	 */
	fprintf(outfp, spaces); /* blank line */
	sprintf(writebuf, "GROUP   1010");
	k = strlen(writebuf);
	strcpy(&writebuf[k], &spaces[k]);
	fprintf(outfp, writebuf);
	fprintf(outfp, spaces); /* blank line */
	/*
	   Header title lines with dates, for example:

	   JPL Planetary Ephemeris DE405/LE405
	   Start Epoch: JED=  2305424.5 1599 DEC 09 00:00:00
	   Final Epoch: JED=  2525008.5 2201 FEB 20 00:00:00
	 */
	sprintf(header->ttl[0],"JPL Planetary Ephemeris DE%03d/LE%03d",
			header->numde, header->numle);
	k = strlen(header->ttl[0]);
	strcpy(&header->ttl[0][k], &spaces[k]);
	ephcom_jd2cal(header->ss[0], idate, 0);
	sprintf(header->ttl[1],"Start Epoch: JED=%11.1f%5d %3s %02d %02d:%02d:%02d",
			header->ss[0], idate[0], month[idate[1]-1], idate[2],
			idate[3], idate[4], idate[5]);
	k = strlen(header->ttl[1]);
	strcpy(&header->ttl[1][k], &spaces[k]);
	ephcom_jd2cal(header->ss[1], idate, 0);
	sprintf(header->ttl[2],"Final Epoch: JED=%11.1f%5d %3s %02d %02d:%02d:%02d",
			header->ss[1], idate[0], month[idate[1]-1], idate[2],
			idate[3], idate[4], idate[5]);
	k = strlen(header->ttl[2]);
	strcpy(&header->ttl[2][k], &spaces[k]);

	/*
	   Don't print trailing blanks at the end of these 3 lines.
	 */
	for (i=0; i<3; i++) {
		// strncpy(writebuf, header->ttl[i], 85);
		// for (j=83; isspace(writebuf[j]) || writebuf[j]=='\0'; j--) writebuf[j]='\0';
		// if (i > 0) writebuf[++j] = ' '; /* To match end space in JPL header */
		strncpy(writebuf, header->ttl[i], 80);
		writebuf[80] = '\r';
		writebuf[81] = '\n';
		writebuf[82] = '\0';
		fprintf(outfp, writebuf);
	}
	/*
	   GROUP 1030: Start and End JD, timestep (in JD) per block.
	 */
	fprintf(outfp, spaces); /* blank line */
	sprintf(writebuf, "GROUP   1030");
	k = strlen(writebuf);
	strcpy(&writebuf[k], &spaces[k]);
	fprintf(outfp, writebuf);
	fprintf(outfp, spaces); /* blank line */

	sprintf(writebuf, "%12.2f%12.2f%11.0f.",
			header->ss[0], header->ss[1], header->ss[2]);
	k = strlen(writebuf);
	strcpy(&writebuf[k], &spaces[k]); /* pad with spaces */
	fprintf(outfp, writebuf);
	/*
	   GROUP 1040: Constant names.
	 */
	fprintf(outfp, spaces); /* blank line */
	sprintf(writebuf, "GROUP   1040");
	k = strlen(writebuf);
	strcpy(&writebuf[k], &spaces[k]);
	fprintf(outfp, writebuf);
	fprintf(outfp, spaces); /* blank line */

	sprintf(writebuf, "%6d", header->ncon);
	k = strlen(writebuf);
	strcpy(&writebuf[k], &spaces[k]);
	fprintf(outfp, writebuf);

	/*
	   Now write the constant names, 10 per line, each 6 characters long
	   preceded by 2 blanks.  Pad names with blanks to make 6 characters.
	 */
	for (i=0; i<header->ncon; i++) {
		fprintf(outfp, "  %-6s", header->cnam[i]);
		if (i % 10 == 9) fprintf(outfp, "\r\n");
	}
	if (i % 10 != 0) {  /* Pad last line with spaces (i is 1 more than above) */
		for ( ; i % 10 != 0; i++) fprintf(outfp, "        ");
		fprintf(outfp, "\r\n");
	}
	/*
	   GROUP 1041: Constant values.
	 */
	fprintf(outfp, spaces); /* blank line */
	sprintf(writebuf, "GROUP   1041");
	k = strlen(writebuf);
	strcpy(&writebuf[k], &spaces[k]);
	fprintf(outfp, writebuf);
	fprintf(outfp, spaces); /* blank line */

	sprintf(writebuf, "%6d", header->nval);
	k = strlen(writebuf);
	strcpy(&writebuf[k], &spaces[k]);
	fprintf(outfp, writebuf);

	if (header->nval != header->ncon) {
		fprintf(stderr,"Error: number of constants and values not equal.\n\n");
		exit(1);
	}
	/*
	   Now read constant values, 3 per line, 26 characters each.
	 */
	for (i = 0; i < header->ncon; i += 3) {
		val1 = header->cval[i];
		val2 = (i+1 < header->ncon) ? header->cval[i+1] : 0.0;
		val3 = (i+2 < header->ncon) ? header->cval[i+2] : 0.0;
		sprintf(writebuf, "%25.17E %25.17E %25.17E   \r\n", val1, val2, val3);
		/* Note that the character holding the sign for each # is left as is. */
		ephcom_doublestrc2f(&writebuf[01]);
		ephcom_doublestrc2f(&writebuf[27]);
		ephcom_doublestrc2f(&writebuf[53]);
		fprintf(outfp, "%s", writebuf);
	}
	/*
	   GROUP 1050: Constant values.
	 */
	fprintf(outfp, spaces); /* blank line */
	sprintf(writebuf, "GROUP   1050");
	k = strlen(writebuf);
	strcpy(&writebuf[k], &spaces[k]);
	fprintf(outfp, writebuf);
	fprintf(outfp, spaces); /* blank line */
	/*
	   If there are no coefficients for an ipt[i][] object (i.e., ipt[i][1]==0),
	   then ipt[i][0] should contain the value of the next available coefficient
	   number rather than 0, as per communication of Myles Standish to Paul Hardy
	   on preferred format of ephemeris headers.

	   If there are no libration coefficients (i.e., lpt[1]==0), then lpt[0]
	   should contain the value of the next available coefficient number rather
	   than 0 as well, as per the same communication from Myles Standish.
	 */
	/* First set j to maximum index into ipt[] that has coefficients */
	j = 0;
	for (i=1; i<12; i++)
		if (header->ipt[i][1] > 0 && header->ipt[i][0] > j)
			j = i;
	/* Now set j to next available index count. */
	if (header->lpt[1] > 0 && header->lpt[0] > j)
		j = header->lpt[1] + header->lpt[1] * header->lpt[2] * 3;
	else
		j = header->ipt[j][0] +
			header->ipt[j][1] * header->ipt[j][2] * (j==11 ? 2 : 3);
	for (i=1; i<12; i++)
		if (header->ipt[i][0] == 0) header->ipt[i][0] = j;
	if (header->lpt[0] == 0) header->lpt[0] = j;
	/*
	   Write ipt array in transposed order (arrays are transposed in FORTRAN
	   from their order in C).
	 */
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 12; j++) {
			fprintf(outfp, "%6d", header->ipt[j][i]);
		}
		fprintf(outfp, "%6d  \r\n", header->lpt[i]);
	}
	/*
	   GROUP 1070: Constant values.
	 */
	fprintf(outfp, spaces); /* blank line */
	sprintf(writebuf, "GROUP   1070");
	k = strlen(writebuf);
	strcpy(&writebuf[k], &spaces[k]);
	fprintf(outfp, writebuf);
	fprintf(outfp, spaces); /* blank line */
	/*
	   Now we're pointing to the first block of coefficient data, after header.
	 */
	return(0);
}


/*
   Write coefficient block information in ASCII format.
 */
ephcom_writeascii_block(
		FILE *outfp,
		struct ephcom_Header *header,
		int blocknum,
		double *datablock) {

	double val1, val2, val3; /* To write text line with 3 double precision words */
	int i, j, k, n;

	char writebuf[EPHCOM_MAXLINE + 1];
	char outhcars[EPHCOM_MAXLINE + 1];
	size_t fwrite(const void *ptr, size_t size, size_t  nmemb, FILE *stream);
	int fputc(int, FILE *);
	double ephcom_indouble(FILE *);
	int ephcom_doublestrc2f(char *); /* Convert C formatted double to FORTRAN format */

	/*
	   Write first line in block, which is block number and ncoeff.
	   Note that lines in the data block files end in "\r\n", while
	   lines in the header files just end in "\n".
	 */
	fprintf(outfp, "%6d%6d", blocknum + 1, header->ncoeff);
	for (i=0; i<68; i++) fputc(' ', outfp);
	fprintf(outfp, "\r\n");
	/*
	   Now write the data, 3 coefficients per line, 26 characters each.
	   Convert format to match what appears in JPL Ephemeris ASCII files.
	 */
	for (i = 0; i < header->ncoeff; i += 3) {
		val1 = datablock[i];
		val2 = (i+1) < header->ncoeff ? datablock[i+1] : 0.0;
		val3 = (i+2) < header->ncoeff ? datablock[i+2] : 0.0;
		/*
		   Write values, 3 coefficients per line, pad lines with 0.0E+00
		 */
		sprintf(writebuf,
				"%25.17E %25.17E %25.17E   \r\n",
				val1, val2, val3);
		// printf("[%d]%s==>", strlen(writebuf), writebuf);
		/*
		   Now re-format numbers the way the JPL header file writes them:
		   all with a leading "0.", so the exponent is one greater.
		   Note that here we start at strlen(writebuf)-6, but in the section
		   that handles header data coefficients we start at strlen(writebuf)-5.
		   This is because the data blocks end ASCII lines with "\r\n" whereas
		   the header ends ASCII lines with just "\n".
		 */
		ephcom_doublestrc2f(&writebuf[01]);  /* Reformat first number */
		ephcom_doublestrc2f(&writebuf[27]);  /* Reformat second number */
		ephcom_doublestrc2f(&writebuf[53]);  /* Reformat third number */
		// printf("[%d]%s", strlen(writebuf), writebuf);
		fprintf(outfp, "%s", writebuf);
	}
	return(0);
}

/*
   Write a JPL Ephemeris header in binary format.
 */
int ephcom_writebinary_header(FILE *outfp, struct ephcom_Header *header) {

	char readbuf[EPHCOM_MAXLINE + 1];
	char group[13];  /* To hold "GROUP" header line */
	int blockout; /* number of bytes we've written to current block/rec in file */
	int blockbytes; /* number of bytes in a block, equals 8 * ncoeff */

	double val1, val2, val3; /* To read text line with 3 double precision words */
	int i, j, k, n;
	int idate[6];
	char *month[12] = {"JAN","FEB","MAR","APR","MAY","JUN",
		"JUL","AUG","SEP","OCT","NOV","DEC"};

	char outhcars[EPHCOM_MAXLINE + 1];
	char *ephcom_nxtgrp(char *group, char *expected, FILE *infile);
	size_t fwrite(const void *ptr, size_t size, size_t  nmemb, FILE *stream);

	/*
	   Point to beginning of output file.
	 */
	rewind(outfp);
	/*
	   First header line: KSIZE= # NCOEFF= #
	 */
	blockbytes = sizeof(double) * header->ncoeff;
	/*
	   Start writing output ephemeris, beginning with header.
	 */
	/*
	   Header title lines with dates, for example:

	   JPL Planetary Ephemeris DE405/LE405
	   Start Epoch: JED=  2305424.5 1599 DEC 09 00:00:00
	   Final Epoch: JED=  2525008.5 2201 FEB 20 00:00:00
	 */
	sprintf(header->ttl[0],"JPL Planetary Ephemeris DE%03d/LE%03d",
			header->numde, header->numle);
	for (i=strlen(header->ttl[0]); i<84; i++) header->ttl[1][i] = ' ';
	ephcom_jd2cal(header->ss[0], idate, 0);
	sprintf(header->ttl[1],"Start Epoch: JED=%11.1f%5d %3s %02d %02d:%02d:%02d",
			header->ss[0], idate[0], month[idate[1]-1], idate[2],
			idate[3], idate[4], idate[5]);
	for (i=strlen(header->ttl[1]); i<84; i++) header->ttl[1][i] = ' ';
	ephcom_jd2cal(header->ss[1], idate, 0);
	sprintf(header->ttl[2],"Final Epoch: JED=%11.1f%5d %3s %02d %02d:%02d:%02d",
			header->ss[1], idate[0], month[idate[1]-1], idate[2],
			idate[3], idate[4], idate[5]);
	for (i=strlen(header->ttl[2]); i<84; i++) header->ttl[2][i] = ' ';
	header->ttl[0][84] = header->ttl[1][84] = header->ttl[2][84] = '\0';

	/*
	   ephcom_Header title lines.

	   Write the three title lines to the output file, padded with blanks,
	   84 characters long (84 is the first even multiple of 6 that is > 80,
	   so the JPL software uses that value because it reads in Fortran 'A6'
	   character strings.
	 */
	fprintf(outfp, "%-84s%-84s%-84s", header->ttl[0], header->ttl[1], header->ttl[2]);
	blockout = 3*84;  /* Just wrote 3 84-byte strings to start output file */
	/*
	   Now output 400 cnam entries to the output file.
	 */
	for (i = 0; i < header->ncon; i++) {
		fprintf(outfp, "%-6s", header->cnam[i]);
		blockout += 6;
	}
	for ( ; i < 400; i++) {
		fprintf(outfp, "      ");  /* Round out to 400 entries, all blank at end */
		blockout += 6;
	}
	/*
	   Binary values: Make sure bytes are in big-endian (network) order for file.
	 */
	for (i=0; i<3; i++) {
		ephcom_outdouble(outfp, header->ss[i]);  /* Write net-order bytes from double precision */
		blockout += 8;
	}
	ephcom_outint(outfp, header->ncon);
	blockout += 4;
	ephcom_outdouble(outfp, header->au);
	blockout += 8;
	ephcom_outdouble(outfp, header->emrat);
	blockout += 8;
	/*
	   If there are no coefficients for an ipt[i][] object (i.e., ipt[i][1]==0),
	   then ipt[i][0] should contain the value of the next available coefficient
	   number rather than 0, as per communication of Myles Standish to Paul Hardy
	   on preferred format of ephemeris headers.

	   If there are no libration coefficients (i.e., lpt[1]==0), then lpt[0]
	   should contain the value of the next available coefficient number rather
	   than 0 as well, as per the same communication from Myles Standish.
	 */
	/* First set j to maximum index into ipt[] that has coefficients */
	j = 0;
	for (i=1; i<12; i++)
		if (header->ipt[i][1] > 0 && header->ipt[i][0] > j)
			j = i;
	/* Now set j to next available index count. */
	if (header->lpt[1] > 0 && header->lpt[0] > j)
		j = header->lpt[1] + header->lpt[1] * header->lpt[2] * 3;
	else
		j = header->ipt[j][0] +
			header->ipt[j][1] * header->ipt[j][2] * (j==11 ? 2 : 3);
	for (i=1; i<12; i++)
		if (header->ipt[i][0] == 0) header->ipt[i][0] = j;
	if (header->lpt[0] == 0) header->lpt[0] = j;

	for (j = 0; j < 12; j++) {
		for (i = 0; i < 3; i++) {
			ephcom_outint(outfp, header->ipt[j][i]);
			blockout += 4;
		}
	}
	ephcom_outint(outfp, header->numde);
	blockout += 4;
	for (i = 0; i < 3; i++) {
		ephcom_outint(outfp, header->lpt[i]);
		blockout += 4;
	}
	/*
	   Now pad the end of the first record with null bytes.  Note: the
	   JPL Fortran software just skips to next record at this point.
	 */
	for (i = blockout; i < blockbytes; i++) {
		fputc('\0', outfp);
	}
	/*
	   End of first block.  Now set blockout to 0 and start with next block.
	 */
	blockout = 0;
	for (i=0; i<header->ncon; i++) {
		ephcom_outdouble(outfp, header->cval[i]);
		blockout += 8;
	}
	/*
	   Pad with double-precision zeroes for rest of array.
	 */
	for ( ; i < 400; i++) {
		ephcom_outdouble(outfp, (double)0.0);
		blockout += 8;
	}
	/*
	   Pad with nulls for rest of block.
	 */
	for (i = blockout; i < blockbytes; i++) {
		fputc('\0', outfp);
	}
	/*
	   Finished normally.
	 */
	return(0);
}

/*
   Write a block of data coefficients in JPL binary file format.
 */
ephcom_writebinary_block(
		FILE *outfp,
		struct ephcom_Header *header,
		int blocknum,
		double *datablock) {

	int i;
	int filebyte;
	int filepos;

	/*
	   Find out where we need to point in the binary file.
	 */
	filebyte = (blocknum + 2) * header->ncoeff * 8; /* 8 bytes per coefficient */
	/*
	   If the file isn't that large, pad it with null bytes
	 */
	fseek(outfp, 0L, SEEK_END);
	filepos = ftell(outfp);
	if (filepos < filebyte) {
		for (i=0; i < (filebyte - filepos); i++) {
			fputc('\0', outfp);
		}
	}
	/*
	   Now go to position where we want to start writing.
	 */
	fseek(outfp, filebyte, SEEK_SET);
	for (i = 0; i < header->ncoeff; i++) {
		ephcom_outdouble(outfp, datablock[i]);
	}

	return(0);
}


/*
   ephcom_parse_block() - Parse a binary block of data.  Warning: verbose!
   Writes parsed output to file pointer outfp.
 */
ephcom_parse_block(
		FILE *outfp,
		struct ephcom_Header *header,
		double *datablock
		) {

	int i0, i1, i2, i3;
	int blockword;
	/*
	   Names of the objects in Chebyshev coefficient arrays.
	 */
	static char *ephcom_coeffname[13] = {
		"Mercury","Venus","EMBary","Mars","Jupiter","Saturn","Uranus","Neptune",
		"Pluto","Moon","Sun","Nutation","Libration"};

	blockword=0;
	fprintf(outfp, "@%04d StartJD\t%12.2f\n", blockword++, datablock[0]);
	fprintf(outfp, "@%04d EndJD\t%12.2f\n", blockword++, datablock[1]);
	for (i0=0; i0<13; i0++) {  /* For all bodies */
		fprintf(outfp, "Body\t%d (%s)\n", i0+1, ephcom_coeffname[i0]);
		for (i1=0; i1<header->ipt[i0][2]; i1++) { /* For all subintervals */
			fprintf(outfp, "  Subinterval %d of %d\n", i1+1, header->ipt[i0][2]);
			for (i2=0; i2 < (i0==11 ? 2 : 3); i2++) {  /* For all coordinates */
				fprintf(outfp, "    %cCoefficients\n", 'X' + i2);
				for (i3=0; i3<header->ipt[i0][1]; i3++) {  /* For all coefficients */
					blockword =     header->ipt[i0][0] +
						i1*header->ipt[i0][1] * (i0==11 ? 2 : 3) +
						i2*header->ipt[i0][1] + i3 - 1;
					fprintf(outfp, "      @%04d [%2d of %2d] %25.17E\n",
							blockword, i3+1, header->ipt[i0][1], datablock[blockword]);
				}
			}
		}
	}
	return(0);
}


/*
   Get the next group header in the file.
   Group header lines are 12 characters long, of the form "GROUP   nnnn".

Parameters:
group - the GROUP header we read
expected - the header we expected
infile - the file pointer to read
 */
char *ephcom_nxtgrp(char *group, char *expected, FILE *infile) {

	char readbuf[EPHCOM_MAXLINE + 1];
	char *fgets(char *, int, FILE *);

	fgets(readbuf, EPHCOM_MAXLINE, infile); /* Blank Line     */
	fgets(readbuf, EPHCOM_MAXLINE, infile); /* "GROUP   dddd\n" */
	strncpy(group, readbuf, 12);
	group[12] = '\0';
	if (strncmp(group, expected, 12) != 0) {
		fprintf(stderr, "Badly formed header; \"%s\" not found.\n\n", expected);
		exit(1);
	}
	fgets(readbuf, EPHCOM_MAXLINE, infile); /* Blank Line    */

	return(NULL);
}

/*
   Print a double precision value to the given file with bytes swapped
   if necessary to match network order (Big Endian).  On Intel 80x86
   the bytes will get swapped, on Motorola or SPARC they won't.
 */
ephcom_outdouble(FILE *outfp, double x) {
	double retval;
	unsigned char ch[8];
	unsigned char *gnulliver64c(unsigned char *);
#ifndef __APPLE__
/* JAZ - this declaration seems to confuse gcc on a mac.*/
	void *memcpy(void *, const void *, size_t);
#endif
	memcpy((void *)ch, (const void *)&x, 8);
	// printf("Got OD: %25.17e => %02x %02x %02x %02x %02x %02x %02x %02x, ",
	//        x, ch[0], ch[1], ch[2], ch[3], ch[4], ch[5], ch[6], ch[7]);
	(void)gnulliver64c(ch);
	// printf("Sending OD: %02x %02x %02x %02x %02x %02x %02x %02x\n",
	//        ch[0], ch[1], ch[2], ch[3], ch[4], ch[5], ch[6], ch[7]);
	fwrite(ch, 1, 8, outfp);
	return(0);
}

/*
   Print an integer (4-byte) value to the given file with bytes swapped
   if necessary to match network order (Big Endian).  On Intel 80x86
   the bytes will get swapped, on Motorola or SPARC they won't.
 */
ephcom_outint(FILE *outfp, unsigned u) {
	unsigned u2;
	unsigned gnulliver32(unsigned);

	u2 = gnulliver32(u);
	fwrite(&u2, 4, 1, outfp);
	return(0);
}

/*
   Read in a double precision value from the given file with bytes swapped
   if necessary to match host order (Little- or DEC- Endian).  On Intel 80x86
   the bytes will get swapped, on Motorola or SPARC they won't.
 */
double ephcom_indouble(FILE *infp) {
	double x;
	static double retval;
	size_t fread(void *ptr,  size_t size, size_t nmemb, FILE *stream);
	unsigned char ch[8];
	unsigned char *gnulliver64c(unsigned char *);
#ifndef __APPLE__
/* JAZ - this declaration seems to confuse gcc on a mac.*/
	void *memcpy(void *, const void *, size_t);
#endif

	/*
	   Handle as character string until bytes are in correct order,
	   then copy to double once they are.
	 */
	fread(ch, 1, 8, infp);
	// printf("Got ID: %02x %02x %02x %02x %02x %02x %02x %02x, ",
	//        ch[0], ch[1], ch[2], ch[3], ch[4], ch[5], ch[6], ch[7]);
	(void)gnulliver64c(ch);
	memcpy((void *)&retval, (const void *)ch, (size_t)8);
	// printf("Sending ID: %02x %02x %02x %02x %02x %02x %02x %02x [%25.18E]\n",
	//        ch[0], ch[1], ch[2], ch[3], ch[4], ch[5], ch[6], ch[7], retval);
	return(retval);
}

/*
   Read in an integer (4--byte) value to the given file with bytes swapped
   if necessary to match host order (Little- or DEC- Endian).  On Intel 80x86
   the bytes will get swapped, on Motorola or SPARC they won't.
 */
int ephcom_inint(FILE *infp) {
	unsigned u;
	static int retval;
	unsigned gnulliver32(unsigned);

	fread(&u, 4, 1, infp);
	retval = (int)gnulliver32(u);
	return(retval);
}


/*
   ephcom_doublstrc2f() - function to convert a string with a double precision
   value written in C to a double precision value that
   FORTRAN creates.  Conversion happens in place.
 */

ephcom_doublestrc2f(char *buf) {

	int i, j, istart, istop, exp, edigits;
	double x;

	for (istart=0; isspace(buf[istart]); istart++);
	x = atof(&buf[istart]);
	// printf("x = %E\n", x);
	for (istop=istart; toupper(buf[istop]) != 'E'; istop++);
	exp = atoi(&buf[istop+1]);
	exp++;
	if (exp < 0) {
		buf[istop+2] = '-';
		exp = -exp;
	}
	else {
		buf[istop+2] = '+';
	}
	if (x == 0.0) exp=0;
	if (exp < 100) edigits = 2;
	else if (exp < 1000) edigits = 3;
	else edigits = 4;
	// printf("istart=%d, istop=%d, exp=%d, edigits=%d\n",
	//        istart, istop, exp, edigits);

	while (edigits > 0) {
		buf[istop + edigits + 2] = exp % 10 + '0';
		exp /= 10;
		edigits--;
	}

	buf[istop+1] = 'D';

	while (istop > istart && buf[istop-1] != '.') {
		buf[istop] = buf[istop-1];
		istop--;
	}

	buf[istop] = buf[istop-2];  /* buf[istop-1] == '.' */
	buf[istop-2] = '0';         /* leading zero */

	return(0);
}


/*
   Planetary Ephemeris.  Takes coordinates already calculated in
   coords structure an converts to vectors and vector dot in testr[].
   Bodies start at 1 for Mercury, to match the JPL PLEPH() numbering.
   Values for ntarg and ncntr correspond to locations ntarg-1 and
   ncntr-1 in coords->pv[].
 */
ephcom_pleph(struct ephcom_Coords *coords, int ntarg, int ncntr, double *r) {

	int i,j;

	if (ntarg != 14 && ncntr != 14) { /* If not nutation, handle normally */
		if (ntarg == 15 || ncntr == 15) { /* Libration */
			for (i=0; i<6; i++) r[i] = coords->pv[14][i];
			//    for (i=0; i<6; i++) printf("\nlibration[%d] = %g\n", i, r[i]);
		}
		else {
			for (i=0; i<6; i++) {
				r[i] = coords->pv[ntarg-1][i] - coords->pv[ncntr-1][i];
				//       printf("\n%g = %g - %g\n",
				//              r[i], coords->pv[ntarg-1][i], coords->pv[ncntr-1][i]);
			}
		}
	}
	else { /* Nutation */
		r[0] = coords->pv[13][0];
		r[1] = coords->pv[13][1];
		r[2] = coords->pv[13][2];
		r[3] = coords->pv[13][3];
		r[4] = 0.0;
		r[5] = 0.0;
		// for (i=0; i<4; i++) printf("\nnutation[%d] = %g\n", i, r[i]);
	}

	return(0);
}


/*
   ephcom_get_coords() - Interpolate positions and velocities at given time.
 */
ephcom_get_coords(FILE *infp,
		struct ephcom_Header *header,
		struct ephcom_Coords *coords,
		double *datablock) {

	double et2[2];    /* Ephemeris time, as coarse (whole) and fine time  in JD */
	double totaltime; /* Sum of whole and fractional JD */
	double filetime;  /* JDs since start of ephemeris file */
	double blocktime; /* JDs since start of data block */
	double subtime;   /* JDs since start of subinterval in block */

	int i, j, k;
	int blocknum;
	int nsub; /* Number of subintervals in data block for this body */
	int subinterval; /* Number of subinterval for this body */
	int dataoffset; /* Offset in datablock for current body and subinterval */
	double subspan; /* Span of one subinterval in days */
	double chebytime; /* Normalized Chebyshev time, in interval [-1,1]. */
	int ncoords; /* Number of coordinates for position and velocity */
	int ncf; /* Number of Chebyshev coefficients per coordinate */
	int retval; /* Return value */

	retval = 0; /* Assume normal return */
	/*
	   Split time JD into whole JDs (et2[0]) and fractional JD (et2[1]).
	 */
	totaltime = coords->et2[0] + coords->et2[1];
	if (totaltime < header->ss[0] || totaltime > header->ss[1]) {
		// fprintf(stderr,"Time is outside ephemeris range.\n");
		retval = -1;
	}
	else {
		et2[0] = (int)totaltime;
		et2[1] = (coords->et2[0] - et2[0]) + coords->et2[1];
		filetime = totaltime - header->ss[0]; /* Days from start of file */
		blocknum = (int)(filetime / header->ss[2]); /* Data block in file, 0.. */
		/*
		   Read the data block that contains coefficients for desired date
		 */
		ephcom_readbinary_block(infp, header, blocknum, datablock);
		// printf("%12.2f : Block %d, %12.2f to %12.2f\n",
		//        totaltime, blocknum, datablock[0], datablock[1]);
		// ephcom_parse_block(stdout, header, datablock);
		/*
		   Now step through the bodies and interpolate positions and velocities.
		 */
		blocktime = totaltime - datablock[0]; /* Days from block start */
		for (i=0; i<13; i++) {
			//    printf("Calculating Body %d\n", i+1);
			subspan = header->ss[2] / header->ipt[i][2]; /* Days/subinterval */
			subinterval = (int)((totaltime - datablock[0]) / subspan);
			if (coords->seconds) subspan *= 86400.0; /* for km/sec, 86400 sec/day */
			//    printf("subspan = %g\n", subspan);
			ncoords = i == 11 ? 2 : 3; /* 2 coords for nutation, else 3 */
			//    printf("Body %d for JD %12.2f subinterval=%d [%d %d %d]\n",
			//           i, totaltime, subinterval,
			//              header->ipt[i][0], header->ipt[i][1], header->ipt[i][2]);
			dataoffset = header->ipt[i][0] - 1 +
				ncoords * header->ipt[i][1] * subinterval;
			//    printf("Body %d, Subinterval %d, Dataoffset = %d\n",
			//           i+1, subinterval, dataoffset);
			subtime = blocktime - subinterval * header->ss[2] / header->ipt[i][2];
			/*
			   Divide days in this subblock by total days in subblock
			   to get interval [0,1].  The right part of the expression
			   will evaluate to a whole number: subinterval lengths are
			   all integer multiples of days in a block (all powers of 2).
			 */
			chebytime = subtime / subspan;
			chebytime = (chebytime + chebytime) - 1.0;
			if (chebytime < -1.0 || chebytime > 1.0) {
				fprintf(stderr,"Chebyshev time is beyond [-1,1] interval.\n");
				fprintf(stderr,
						"filetime=%f, blocktime=%f, subtime=%f, chebytime=%f\n",
						filetime, blocktime, subtime, chebytime);
			}
			else {
				/*
				   Everything is as expected.  Interpolate coefficients.
				 */
				ephcom_cheby(header->maxcheby, chebytime, subspan,
						&datablock[dataoffset],
						ncoords, header->ipt[i][1], coords->pv[i]);
			}
		}
		// /*
		//    Set any user-defined coordinates to zero.
		// */
		// for (i = 16; i < EPHCOM_NUMOBJECTS; i++)
		//    coords->pv[i][0] = coords->pv[i][1] = coords->pv[i][1] =
		//       coords->pv[i][1] = coords->pv[i][1] = coords->pv[i][1] =  0.0;
		/*
		   With interpolations complete, calculate Earth from EMBary and
		   Sun from SSBary.  Preserve other coordinates.
		 */
		for (j=0; j<6; j++) {
			coords->pv[15][j] = coords->pv[ 9][j]; /* Save original lunar coords */
			coords->pv[14][j] = coords->pv[12][j]; /* Librations if on file */
			coords->pv[13][j] = coords->pv[11][j]; /* Nutations if on file */
			/*
			   Calculate Earth and Moon from EMBary and geocentric Moon.
			 */
			coords->pv[12][j] = coords->pv[2][j]; /* Move EMBary from Earth spot */
			coords->pv[2][j] -= coords->pv[9][j] / (1.0 + header->emrat); /* Earth */
			coords->pv[9][j] += coords->pv[2][j]; /* Moon (change geo->SS-centric) */
		}
		/*
		   If we want heliocentric coordinates (not Solar System Barycenter
		   coordinates), subtract Sun's position from all coordinates.
		 */
		if (!coords->bary) {
			for (i=0; i<13; i++) {
				if (i == coords->center) i++; /* Skip over center object */
				if (i < 13) {
					for (j=0; j<6; j++)
						coords->pv[i][j] -= coords->pv[coords->center][j];
				}
			}
			/*
			   Set new center's positions and velocities to all zeroes
			 */
			coords->pv[coords->center][0] = coords->pv[coords->center][1] =
				coords->pv[coords->center][2] = coords->pv[coords->center][3] =
				coords->pv[coords->center][4] = coords->pv[coords->center][5] = 0.0;
		}
		if (!coords->km) { /* Calculate AU, not kilometers */
			for (i=0; i<15; i++ ) {
				if (i == 13) i = 15; /* Skip over nutations and librations */
				for (j=0; j<6; j++)
					coords->pv[i][j] /= header->au;
			}
		}
	}

	return(retval);
}


/*
   ephcom_cheby() - interpolate at a point using Chebyshev coefficients
 */
inline ephcom_cheby(
		int maxcoeffs, /* Maximum number of Chebyshev components possible */
		double x,      /* Value of x over [-1,1] for Chebyshev interpolation */
		double span,   /* Span in time of subinterval, for velocity */
		double *y,     /* Chebyshev coefficients */
		int ncoords,   /* Total number of coordinates to interpolate */
		int ncoeffs,   /* Number of Chebyshev coefficients per coordinate */
		double *pv     /* Array to hold position in 1st half, velocity in 2nd */
		) {

	int i, j;
	static double twox;
	static double *pc, *vc; /* Position and velocity polynomial coefficients. */
	static double lastx=2.0; /* x from last call; initialize to impossible value */
	static init=1; /* Need to initialize pc[] and vc[] */

	/*
	   Allocate position and velocity Chebyshev coefficients.
	 */
	if (init) {
		pc = (double *)malloc(maxcoeffs * sizeof(double));
		vc = (double *)malloc(maxcoeffs * sizeof(double));
		init = 0;
	}
	/*
	   This need only be called once for each Julian Date,
	   saving a lot of time initializing polynomial coefficients.
	 */
	if (lastx != x) {
		lastx = x;
		twox = x + x;   /* For Chebyshev recursion */
		/*
		   Initialize position polynomial coefficients
		 */
		pc[0] = 1.0;    /* Chebyshev T[0](x) = 1 */
		pc[1] = x;      /* Chebyshev T[1](x) = x */
		// printf("pc[%d] = %20.13e\n", 0, pc[0]);
		// printf("pc[%d] = %20.13e\n", 1, pc[1]);
		for (i=2; i<maxcoeffs; i++) {
			pc[i] = twox * pc[i-1] - pc[i-2];
			//    printf("pc[%d] = %20.13e\n", i, pc[i]);
			/*
			   Resolve bug with gcc generating -0.0 (also makes
			   the smallest represented number equal to zero).
			 */
			if (pc[i]*pc[i] == 0.0) {
				pc[i] = 0.0;
				//       printf("pc[%d] = %20.13e\n", i, pc[i]);
			}
		}
		/*
		   Initialize derivative polynomial coefficients
		 */
		vc[0] = 0.0;          /* d(1)/dx        = 0  */
		vc[1] = 1.0;          /* d(x)/dx        = 1  */
		vc[2] = twox + twox;  /* d(2x^2 - 1)/dx = 4x */
		for (i=3; i<maxcoeffs; i++) {
			vc[i] = twox * vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2];
			//    printf("vc[%2d] = 2(%g) * %g + 2 * %g - %g = %g\n",
			//           i, x, vc[i-1], pc[i-1], vc[i-2], vc[i]);
		}
	}
	/*
	   Interpolate to get position for each component
	 */
	for (i=0; i<ncoords; i++) { /* Once each for x, y, and z */
		pv[i] = 0.0;
		for (j=ncoeffs-1; j >= 0; j--) {
			pv[i] += pc[j] * y[i*ncoeffs + j];
			//    printf("[%2d;%4d]pv[%2d] += %17.10e * %17.10e = %17.10e\n",
			//           j, i*ncoeffs+j, i, pc[j], y[i*ncoeffs+j], pv[i]);
		}
	}
	/*
	   Interpolate velocity (first derivative)
	 */
	for (i=0; i<ncoords; i++) {
		pv[ncoords + i] = 0.0;
		// printf("pv[%2d] = 0.0\n", ncoords+i);
		for (j=ncoeffs-1; j >= 0; j--) {
			pv[ncoords + i] += vc[j] * y[i*ncoeffs + j];
			//    printf("[%2d;%4d]pv[%d] += %g * %17.10e = %17.10e\n",
			//           j, i*ncoeffs+j, ncoords+i, vc[j], y[i*ncoeffs+j], pv[ncoords+i]);
		}
		// printf("[%2d]pv[%d] *= 2.0 / %17.10e = ", i, ncoords+i, span);
		pv[ncoords + i] *= 2.0 / span;
		// printf("%17.10e\n", pv[ncoords+i]);
	}
	return(0);
}

/*
   ephcom_jd2cal() - convert Julian Day to calendar date and time.

tjd: double precision Julian Day
idate: integer year, month, day, hour, minute, second of tjd
calendar_type: -1=Julian; 0=Automatic; 1=Gregorian

If automatic, use Julian calendar for dates before 15 October 1582.

From pp. 604, 606 in the Explanatory Supplement to the Astronomical Almanac.
 */
ephcom_jd2cal(double tjd, int idate[6], int calendar_type) {

	int ihour, imin, isec;
	int j;
	/* From Explanatory Supplement to Astronomical Almanac, pp. 604, 606 */
	int I, J, K, L, N, D, M, Y;

	tjd += 0.5 + 0.5/86400.0; /* Round to nearest second */
	j = tjd;  /* Integer Julian Day */
	tjd = (tjd - j) * 24.0;
	ihour = tjd;
	tjd = (tjd - ihour) * 60.0;
	imin = tjd;
	tjd = (tjd - imin) * 60.0;
	isec = tjd;
	/*
	   Julian calendar.  Explanatory Supplement to Astronomical Alamanac, p. 606.
	   If automatic, use Julian calendar for dates before 15 October 1582.
	 */
	if (calendar_type == -1 || (calendar_type == 0 && j <= 2299160)) {
		J = j + 1402;
		K = (J - 1) / 1461;
		L = J - 1461 * K;
		N = (L - 1) / 365 - L / 1461;
		I = L - 365 * N + 30;
		J = (80 * I) / 2447;
		D = I - (2447 * J) / 80;
		I = J / 11;
		M = J + 2 - 12 * I;
		Y = 4 * K + N + I - 4716;
	}
	/*
	   Gregorian calendar.
	 */
	else { /* Explanatory Supplement to Astronomical Almanac, p. 604 */
		L = j + 68569;
		N = (4 * L) / 146097;
		L = L - (146097 * N + 3) / 4;
		I = (4000 * (L + 1)) / 1461001;
		L = L - (1461 * I) / 4 + 31;
		J = (80 * L) / 2447;
		D = L - (2447 * J) / 80;
		L = J / 11;
		M = J + 2 - 12 * L;
		Y = 100 * (N - 49) + I + L;
	}

	idate[0] = Y;
	idate[1] = M;
	idate[2] = D;
	idate[3] = ihour;
	idate[4] = imin;
	idate[5] = isec;

	return(0);
}

/*
   ephcom_cal2jd() - convert calendar date and time to JD.

idate: integer year, month, day, hour, minute, second
calendar_type: -1=Julian; 0=Automatic; 1=Gregorian
return value: double precision Julian Day of idate[]

From pp. 604, 606 in the Explanatory Supplement to the Astronomical Almanac.
 */
double ephcom_cal2jd(int idate[6], int calendar_type) {

	double tjd;
	int jd;

	/*
	   Convert hours, minutes, seconds to fractional JD.
	 */
	tjd = (idate[3] + (idate[4] + idate[5] / 60.0) / 60.0) / 24.0 - 0.5;
	/*
	   Julian calendar.  Explanatory Supplement to Astronomical Alamanac, p. 606.
	   If automatic, use Julian calendar for dates before 15 October 1582.
	 */
	if (calendar_type == -1 ||
			(calendar_type == 0 && 
			 (idate[0] < 1582 ||                         /* Before 1582 */
			  (idate[0] == 1582 &&
			   (idate[1] < 10 ||                         /* Before October 1582 */
				(idate[1] == 10 && idate[2] < 15)))))) { /* Before 15 October 1582 */
		jd = 367 * idate[0] -
			(7 * (idate[0] + 5001 + (idate[1] - 9) / 7)) / 4 +
			(275 * idate[1]) / 9 +
			idate[2] + 1729777;
	}
	/*
	   Gregorian calendar.
	 */
	else { /* Explanatory Supplement to Astronomical Almanac, p. 604 */
		jd = (1461 * (idate[0] + 4800 + (idate[1] - 14) / 12)) / 4 +
			(367 * (idate[1] - 2 - 12 * ((idate[1] - 14) / 12))) / 12 -
			(3 * ((idate[0] + 4900 + (idate[1] - 14) / 12) / 100)) / 4 +
			idate[2] - 32075;
	}
	/*
	   Return value is whole JD number plus fractional JD number.
	 */
	tjd += (double)jd;

	return(tjd);
}
