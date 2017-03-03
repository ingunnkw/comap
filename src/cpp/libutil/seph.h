#ifndef SEPHH
#define SEPHH

#include <vector>
#include <string>
#include <sgrid.h>
extern "C" {
#include <ephcom.h>
}

namespace skn {

extern const double mdiff;
extern const Grid<double,2> equ2gal;

struct Ephfile
{
	Ephfile(const std::string & s);
	void assert_in_range(double jd);
	FILE * db;
	ephcom_Header header;
	std::vector<double> data;
};

int lookup_object(const std::string & name);
Vector<double,3> get_pos(int obj, double jd, Ephfile & eph, int other = EPHCOM_EARTH);
Vector<double,3> get_pos_rel(int obj, double jd, Ephfile & eph);

Grid<double,2> euler_zyz(double phi, double theta, double psi);
Vector<double,3> operator*(const Grid<double,2> & m, const Vector<double,3> & v);
// Polar coordinates are in (lon,lat,r) convention in radians, not in
// the normal mathematical one
Vector<double,3> rec2pol(const Vector<double,3> & rec);
Vector<double,3> pol2rec(const Vector<double,3> & pol);

struct Object {std::string name; Vector<double,2> pos; double size; bool fixed; int priority;};
std::vector<Object> read_objs(const std::string & filename);

bool match(const std::string & a, const std::string & b, int n);
double timediff(const Vector<double,3> & d);

std::string identify(Ephfile & eph, const std::vector<Object> & objs, double mjd, double lon, double lat, double * dist = 0);

}

#endif
