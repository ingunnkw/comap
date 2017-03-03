#include <scolorizer.h>
#include <string>
#include <vector>
#include <svector.h>
#include <stoken.h>
#include <stypes.h>
#include <serror.h>
#include <cstdio>

namespace skn {

Colorizer::Colorizer() {}
Colorizer::Colorizer(const std::string & s)
{
	std::vector<std::string> tokens = tokenize(s, ", ");
	for(int i = 0; i < tokens.size(); i++)
	{
		std::vector<std::string> parts = tokenize(tokens[i], ": #");
		if(parts.size() != 2)
			serror("Syntax serror in colorizer specification '%s'!", tokens[i].c_str());
		double p; Pixel col(0,0,0);
		if(std::sscanf(parts[0].c_str(), "%lf", &p) < 1)
			serror("Error parsing colorization location in '%s'!", parts[0].c_str());
		col = parse_color(parts[1].c_str());
		if(std::sscanf(parts[1].c_str(), "%02hhx%02hhx%02hhx", &col(0), &col(1), &col(2)) < 3)
			serror("Error parsong colorization color in '%s'!", parts[1].c_str());
		pos.push_back(p);
		cols.push_back(col);
	}
}

Colorizer & Colorizer::push_back(double val, const Pixel & col)
{
	pos.push_back(val);
	cols.push_back(col);
	return *this;
}

Pixel Colorizer::operator()(double d) const
{
	if(pos.empty()) return Pixel();
	if(pos.size() == 1) return cols[0];
	if(d > pos.back()) return cols.back();
	if(d < pos.front()) return cols.front();
	int b; for(b = 1; b < pos.size() && pos[b] < d; b++);
	int a = b-1;
	double f = (d-pos[a])/(pos[b]-pos[a]), g = (pos[b]-d)/(pos[b]-pos[a]);
	Pixel res;
	for(int i = 0; i < size(res); i++)
		res[i] = byte(g*cols[a](i) + f*cols[b](i) + 0.5);
	return res;
}

Colorizer Colorizer::wmap("0:000080,0.15:0000ff,0.4:00ffff,0.7:ffff00,0.9:ff5500,1:800000");

Pixel parse_color(const char * str)
{
	Pixel col;
	if(std::sscanf(str, "%02hhx%02hhx%02hhx", &col(0), &col(1), &col(2)) < 3)
		serror("Error parsing colorization color in '%s'!", str);
	return col;
}

}
