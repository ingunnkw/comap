#ifndef STYPEBB
#define STYPEBB

#include <stdint.h>
#include <svector.h>
#include <sgrid.h>

namespace skn
{
	typedef uint8_t byte;
	typedef Vector<byte,3> Pixel;
	typedef Grid<Pixel,2> Image;
}

#endif
