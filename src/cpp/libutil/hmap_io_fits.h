#ifndef HMAPIOFITSH
#define HMAPIOFITSH

/* This header provides read_hmap_fits and write_hmap_fits */

#include <hmap.h>
#include <fitshandle.h>
#include <fits_types.h>
#include <lsconstants.h>

namespace skn
{

std::string mkname_safe(const std::string & s, int d)
{
	if(s.size() == 0) return s;
	for(int n = s.size()+10;; n *= 2)
	{
		std::vector<char> buf(n);
		if(snprintf(&buf[0], n, s.c_str(), d) < n) return std::string(&buf[0]);
	}
}

bool file_exists(const std::string & s) {
	FILE * f = fopen(s.c_str(), "r");
	bool exists = f;
	if(f) fclose(f);
	return f;
}

template<class T>
HMap<T> read_hmap_fits_single(const char * fname, int hdu = 2)
{
	using namespace std;
	//fitshandle f(fname);
	fitshandle f;
	f.open(fname);
	f.goto_hdu(hdu);

	string order;
	f.get_key("ORDERING", order);
	int ordering = order == "RING" ? 1 : 2;

	int nside; f.get_key("NSIDE", nside);
	int npix = 12*nside*nside;
	string indexschm = "IMPLICIT";
	if(f.key_present("INDXSCHM")) f.get_key("INDXSCHM", indexschm);

	Grid<int,1> pixels;
	int off = 0, n, fullsky = 1;
	if(indexschm == "IMPLICIT") n = npix;
	else
	{
		pixels.resize(f.nelems(1));
		f.read_column_raw(1, &pixels[0], f.nelems(1));
		off++;
		n = size(pixels);
		fullsky = 0;
	}
	Grid<T,3> data(1, f.ncols()-off, n);
	for(int i = 0; i < f.ncols()-off; i++)
		f.read_column_raw(1+i+off, &data(0,i,0), data.size(2));
	if(fullsky) return HMap<T>(nside, ordering, data);
	else return HMap<T>(nside, ordering, pixels, data);
}

template <class T>
HMap<T> read_hmap_fits(const std::string & fname, int hdu = 2)
{
	using namespace std;
	vector<HMap<T> > imaps;
	string prev;
	bool series = fname.find('%') != string::npos;
	for(int m = 0;;m++)
	{
		string iname = mkname_safe(fname, m);
		if(iname == prev) break;
		bool exists = file_exists(iname);
		if(!exists) {
			if(series) break;
			else serror("%s does not exist!", iname.c_str());
		}
		imaps.push_back(read_hmap_fits_single<T>(iname.c_str(),hdu));
		prev = iname;
	}
	return concatenate(imaps);
}

template <class T>
void write_hmap_fits(const std::string & fname, const HMap<T> & map, int type)
{
	// Output all maps, or just the first if no % i present
	using namespace std;
	int mlen = map.size(0), npix = 12*map.get_nside()*map.get_nside();
	if(fname.find('%') == string::npos) mlen = 1;
	for(int m = 0; m < mlen; m++)
	{
		string oname = mkname_safe(fname, m);
		fitshandle f; f.create("!" + oname);
		vector<fitscolumn> cs;
		if(!map.is_fullsky())
			cs.push_back(fitscolumn("INDEX","Pixel index",1,PLANCK_INT32));
		for(int i = 0; i < map.size(1); i++)
		{
			char buf[0x100];
			sprintf(buf, "Signal%d", i+1);
			cs.push_back(fitscolumn(buf, "unknown", 1, (PDT)type));
		}
		f.insert_bintab(cs);
		f.set_key("PIXTYPE", string("HEALPIX"),"HEALPIX pixelisation");
		f.set_key("ORDERING",string(map.get_ordering() == 1 ? "RING":"NESTED"),"Pixel ordering scheme, either RING or NESTED");
		f.set_key("NSIDE", map.get_nside(), "Resolution parameter for HEALPIX");
		f.set_key("FIRSTPIX",0,"First pixel # (0 based)");
		f.set_key("LASTPIX",npix,"Last pixel # (0 based)");
		f.set_key("INDXSCHM",string(map.is_fullsky() ? "IMPLICIT" : "EXPLICIT"),
			"Indexint: IMPLICIT or EXPLICIT");

		int n = map.size(2), off = 1-map.is_fullsky();
		if(!map.is_fullsky())
			f.write_column_raw(1, &map.get_pixels()[0], n);
		for(int i = 0; i < map.size(1); i++)
			f.write_column_raw(1+i+off,&map(m,i,0),n);
	}
}

template <class T>
void write_hmap_fits(const std::string & fname, const HMap<T> & map, const char * dtype = 0)
{
	write_hmap_fits(fname, map, dtype ? dtype2fits(dtype) : fitstype<T>());
}

}

#endif
