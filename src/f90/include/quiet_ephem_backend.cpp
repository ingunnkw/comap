#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
extern "C" {
#include <ephcom.h>
}
using namespace std;
const double pi = 4*atan(1);
const double mdiff = 2400000.5;

struct Ephfile
{
	Ephfile(const char * s)
	{
		db = fopen(s, "r");
		if(!db) {
			fprintf(stderr, "Could not open database file %s!", s);
			exit(1);
		}
		ephcom_readbinary_header(db, &header);
		data.resize(header.ncoeff);
	}
	FILE * db;
	ephcom_Header header;
	vector<double> data;
};

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

vector<double> get_pos(int obj, double jd, Ephfile & eph, int other = EPHCOM_EARTH)
{
	ephcom_Coords coords = set_time(jd);
	ephcom_get_coords(eph.db, &eph.header, &coords, &eph.data[0]);

	// And finally get our coordinates
	vector<double> pos_vel(6);
	ephcom_pleph(&coords, obj, other, &pos_vel[0]);
	pos_vel.resize(3);
	return pos_vel;
}

double timediff(const vector<double> & d) { return sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2])/173.144483; }

vector<double> get_pos_rel(int obj, double jd, Ephfile & eph)
{
	double td = 0, diff = 0, ejd = jd, tol = 1e-10;
    int it = 0;
	vector<double> opos, epos = get_pos(EPHCOM_EARTH, jd, eph, EPHCOM_SSBARY);
	do {
		ejd = jd - td;
		opos = get_pos(obj, ejd, eph, EPHCOM_SSBARY);
		for(int i = 0; i < 3; i++) opos[i] -= epos[i];
		td = timediff(opos);
		diff = abs(jd - (ejd + td));
	} while(it++ < 5 && diff > tol);
	return opos;
}

Ephfile * db = NULL;
string db_path = "/projects/quiet/external_data/unix.405";
// Hexagon path
//string db_path = "/work/tonemru/quiet_data/external_data/unix.405";

extern "C" {
  void ephem_equ_rect_c_(int * id, double * mjd, double * out)
  {
    if(!db) db = new Ephfile(db_path.c_str());
    vector<double> p = get_pos_rel(*id, *mjd+mdiff, *db);
    for(int i = 0; i < 3; i++) out[i] = p[i];
  }
  void ephem_equ_rect_abs_c_(int * id, double * mjd, double * out)
  {
    if(!db) db = new Ephfile(db_path.c_str());
    vector<double> p = get_pos(*id, *mjd+mdiff, *db, EPHCOM_SSBARY);
    for(int i = 0; i < 3; i++) out[i] = p[i];
  }
  void ephem_set_db_(char * path, int n)
  {
    db_path = string(path,path+n);
    if(db) { delete db; db = NULL; }
  }
}
