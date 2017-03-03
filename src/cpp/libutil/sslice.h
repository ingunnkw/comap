#ifndef SSLICEH
#define SSLICEH

#include <vector>
#include <string>
#include <svector.h>

namespace skn {

extern const int Any;
typedef std::vector<Vector<int,3> > Slice;

struct cut {
	cut(int from = Any, int num = Any, int step = Any);
	cut & operator()(int from = Any, int num = Any, int step = Any);
	operator const Slice&() const;
	Slice work;
};

Slice parse_slice(const std::string & spec);
Slice apply_slice(const Slice & old, const std::vector<int> lens);

}

#endif
