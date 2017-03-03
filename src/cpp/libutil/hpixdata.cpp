#include <hpixdata.h>
#include <map>
#include <serror.h>
#include <fitshandle.h>
#include <cstdio>

namespace skn {

std::string trypath(const std::vector<std::string> & paths, const std::string & filename)
{
	FILE * f;
	for(int i = 0; i < paths.size(); i++) {
		std::string s = paths[i] + filename;
		if(f = std::fopen(s.c_str(),"r")) {
			fclose(f);
			return s.c_str();
		}
	}
	return std::string();
}
std::vector<std::string> hpix_paths() {
	std::vector<std::string> res;
	res.push_back("./");
	res.push_back("../");
	res.push_back("data/");
	res.push_back("../data/");
	char * hpix = getenv("HEALPIX");
	if(hpix) {
		res.push_back(std::string(hpix) + "/");
		res.push_back(std::string(hpix) + "/data/");
	}
	char * home = getenv("HOME");
	if(home) res.push_back(std::string(home) + "/local/Healpix_2.10/data/");
	return res;
}
arr<double> get_hpix_array(int nside, char * fmt, double add)
{
	static std::vector<std::string> paths = hpix_paths();
	char filename[0x100];
	sprintf(filename, fmt, nside);
	std::string path = trypath(paths, filename);
	if(path.empty()) skn::serror("Could not locate %s", filename);
	//fitshandle f(path);
	fitshandle f;
	f.open(path);
	f.goto_hdu(2);
	arr<double> a(f.nrows());
	f.read_column(1, a);
	for(int i = 0; i < a.size(); i++) a[i] += add;
	return a;
}
const arr<double> & get_weights(int nside, int lmax)
{
	typedef std::map<int, arr<double> > Store;
	static Store store;
	Store::iterator i = store.find(nside);
	if(!lmax) lmax = 2*nside;
	if(i != store.end()) return i->second;
	else
	{
		arr<double> tmp = get_hpix_array(nside,"weight_ring_n%05d.fits",1);
		if(tmp.size() >= lmax) return store[nside] = tmp;
		arr<double> res(lmax);
		int i; for(i = 0; i < tmp.size(); i++) res[i] = tmp[i];
		for(; i < lmax; i++) res[i] = 1;
		return store[nside] = res;
	}
}
const arr<double> & get_window(int nside)
{
	typedef std::map<int, arr<double> > Store;
	static Store store;
	Store::iterator i = store.find(nside);
	if(i != store.end()) return i->second;
	else return store[nside] = get_hpix_array(nside,"pixel_window_n%04d.fits",0);
}

}
