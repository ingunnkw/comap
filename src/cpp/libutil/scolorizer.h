#ifndef SCOLORINCG
#define SCOLORINCG

#include <string>
#include <vector>
#include <stypes.h>

namespace skn {

class Colorizer
{
public:
	Colorizer();
	Colorizer(const std::string & s);
	Colorizer & push_back(double val, const Pixel & col);
	Pixel operator()(double d) const;
	static Colorizer wmap;
private:
	std::vector<double> pos;
	std::vector<Pixel> cols;
};

Pixel parse_color(const char * str);

}

#endif
