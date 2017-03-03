#include <seph.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <serror.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>

namespace skn {

using std::cos;
using std::sin;
using std::atan;
using std::atan2;
using std::sqrt;
using std::string;
using std::vector;
using std::abs;
using std::acos;
using std::ifstream;

const double pi = 4*atan(1);
const double mdiff = 2400000.5;

typedef Grid<double,2> Matrix;
typedef Vector<double,3> Vec;

Matrix euler_zyz(double phi, double theta, double psi)
{
	Matrix euler_matrix(3,3);

	double sphi = sin(phi);
	double cphi = cos(phi);

	double sth  = sin(theta);
	double cth  = cos(theta);

	double spsi = sin(psi);
	double cpsi = cos(psi);

	euler_matrix(0,0) = -sphi * spsi + cth * cphi * cpsi;
	euler_matrix(0,1) = -sphi * cpsi - cth * cphi * spsi;
	euler_matrix(0,2) =                sth * cphi;
	euler_matrix(1,0) =  cphi * spsi + cth * sphi * cpsi;
	euler_matrix(1,1) =  cphi * cpsi - cth * sphi * spsi;
	euler_matrix(1,2) =                sth * sphi;
	euler_matrix(2,0) =              - sth * cpsi;
	euler_matrix(2,1) =                sth * spsi;
	euler_matrix(2,2) =                cth;
	return euler_matrix;
}

const Matrix equ2gal = euler_zyz(-0.996030630370769, -1.09731754119245, -3.36603324368708);

Ephfile::Ephfile(const string & s)
{
	db = fopen(s.c_str(), "r");
	if(!db) serror("Could not open database file %s!", s.c_str());
	ephcom_readbinary_header(db, &header);
	data.resize(header.ncoeff);
}
void Ephfile::assert_in_range(double jd)
{
	if(jd < header.ss[0] || jd > header.ss[1])
		serror("Specified mjd %.4f is outside data range (%.4f:%.4f)!",
			jd-mdiff, header.ss[0]-mdiff, header.ss[1]-mdiff);
}

ephcom_Coords set_time(double jd)
{
	ephcom_Coords coords;
	coords.km = 0;       // AU
	coords.seconds = 0;  // days
	coords.bary = 1;     // Use barycenter internally
	coords.et2[0] = jd;
	coords.et2[1] = 0;
	return coords;
}

bool match(const string & a, const string & b, int n)
{
	return a.size() <= b.size() && a.size() >= n
		&& a == b.substr(0,a.size());
}

int lookup_object(const string & name)
{
	if(match(name,"sun",2)) return EPHCOM_SUN;
	else if(match(name,"mercury",2)) return EPHCOM_MERCURY;
	else if(match(name,"venus",1)) return EPHCOM_VENUS;
	else if(match(name,"earth",1)) return EPHCOM_EARTH;
	else if(match(name,"moon",2)) return EPHCOM_MOON;
	else if(match(name,"mars",2))  return EPHCOM_MARS;
	else if(match(name,"jupiter",1)) return EPHCOM_JUPITER;
	else if(match(name,"saturn",2)) return EPHCOM_SATURN;
	else if(match(name,"uranus",1)) return EPHCOM_URANUS;
	else if(match(name,"neptune",1)) return EPHCOM_NEPTUNE;
	else if(match(name,"pluto",1)) return EPHCOM_PLUTO;
	else serror("Unknown object %s!", name.c_str());
}

Vec get_pos(int obj, double jd, Ephfile & eph, int other)
{
	ephcom_Coords coords = set_time(jd);
	ephcom_get_coords(eph.db, &eph.header, &coords, &eph.data[0]);

	// And finally get our coordinates
	vector<double> pos_vel(6);
	ephcom_pleph(&coords, obj, other, &pos_vel[0]);
	return Vec(pos_vel[0],pos_vel[1],pos_vel[2]);
}

Vec operator*(const Matrix & m, const Vec & v)
{
	Vec res;
	for(int i = 0; i < size(m,0); i++)
		for(int j = 0; j < size(m,1); j++)
			res[i] += m(i,j)*v[j];
	return res;
}

double timediff(const Vec & d) { return len(d)/173.144483; }

Vec get_pos_rel(int obj, double jd, Ephfile & eph)
{
	double td = 0, diff = 0, ejd = jd, tol = 1e-10;
	int maxit = 10, i = 0;
	Vec opos, epos = get_pos(EPHCOM_EARTH, jd, eph, EPHCOM_SSBARY);
	do {
		ejd = jd - td;
		opos = get_pos(obj, ejd, eph, EPHCOM_SSBARY);
		td = timediff(opos-epos);
		diff = abs(jd - (ejd + td));
	} while(diff > tol && ++i < maxit);
	return opos-epos;
}

// Polar coordinates are in (lon,lat,r) convention in radians, not in
// the normal mathematical one
Vec rec2pol(const Vec & rec)
{
	Vec pol;
	double x = rec[0], y = rec[1], z = rec[2];
	double rp = sqrt(x*x+y*y);
	double r = sqrt(rp*rp+z*z);
	double dec = atan(z/rp), ra = atan2(y,x);
	pol[0] = ra; pol[1] = dec; pol[2] = r;
	return pol;
}

Vec pol2rec(const Vec & pol)
{
	return Vec(cos(pol(1))*cos(pol(0)),cos(pol(1))*sin(pol(0)),sin(pol(1)))*pol(2);
}

struct Match {
	int ind, pri; double dist;
	Match(int ind_, int pri_, double dist_):ind(ind_),pri(pri_),dist(dist_) {}
	bool operator<(const Match & m) const { return pri == m.pri ? dist < m.dist : pri > m.pri; }
};

string identify(Ephfile & eph, const vector<Object> & objs, double mjd, double lon, double lat, double * dist)
{
	double margin = 4*pi/180; // approximate size of focal plane

	vector<Match> matches;
	Vec targ_rect = pol2rec(Vec(lon, lat, 1.0));
	for(int i = 0; i < objs.size(); i++)
	{
		Vec rect;
		if(objs[i].fixed) rect = pol2rec(Vec(objs[i].pos[0],objs[i].pos[1],1));
		else
		{
			int id = lookup_object(objs[i].name);
			rect = equ2gal*get_pos_rel(id, mjd+mdiff, eph);
			rect /= len(rect);
		}
		double dist = acos(dot(targ_rect, rect));
		if(dist < objs[i].size+margin) matches.push_back(Match(i, objs[i].priority, dist));
	}
	sort(matches.begin(), matches.end());
	if(matches.empty()) return "unknown";
	else {
		if(dist) *dist = matches[0].dist;
		return objs[matches[0].ind].name;
	}
}

vector<Object> read_objs(const string & filename)
{
	vector<Object> objs;
	ifstream file(filename.c_str());
	string name, lonstr, latstr;
	double size, res;
	int priority;
	while(file >> name >> lonstr >> latstr >> size >> res >> priority)
	{
		Object obj = { name, Vector<double,2>(), size*pi/180, false, priority };
		if(lonstr != "x" && latstr != "x")
		{
			obj.pos = Vector<double,2>(atof(lonstr.c_str())*pi/180,
				atof(latstr.c_str())*pi/180);
			obj.fixed = true;
		}
		if(priority > 0) objs.push_back(obj);
	}
	return objs;
}

}
