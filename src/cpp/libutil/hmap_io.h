#ifndef HMAPIOH
#define HMAPIOH

#include <hmap_io_hdf.h>
#include <hmap_io_fits.h>

namespace skn {

std::string get_ftype(const std::string & fname, const char * ftype)
{
	using std::string;
	int i = fname.find_last_of(".");
	string filetype = ftype ? ftype : i != string::npos ? fname.substr(i+1) : "";
	return filetype;
}

template <class T>
HMap<T> read_hmap(const std::string & fname, const char * ftype = 0, int hdu = 2)
{
	std::string filetype = get_ftype(fname, ftype);
	if(filetype == "hdf") return read_hmap_hdf<T>(fname);
	else if(filetype == "fits") return read_hmap_fits<T>(fname, hdu);
	else serror("Unrecognized HMap file type: %s", filetype.c_str());
	return HMap<T>();
}

template <class T>
void write_hmap(const std::string & fname, const HMap<T> & map, const char * ftype = 0,
	const char * dtype = 0)
{
	std::string filetype = get_ftype(fname, ftype);
	if(filetype == "hdf") write_hmap_hdf(fname, map, dtype);
	else if(filetype == "fits") write_hmap_fits(fname, map, dtype);
	else serror("Unrecognized HMap file type: %s", filetype.c_str());
}

}

#endif
