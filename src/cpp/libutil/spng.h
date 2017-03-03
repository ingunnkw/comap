#ifndef SPNGINCI
#define SPNGINCI

#include <stypes.h>

namespace skn
{
	// Write a 24-bit rgb png in (x,y) ordering
	void write_png(const char * filename, const Image & img);
	Image read_png(const char * filename);
}

#endif
