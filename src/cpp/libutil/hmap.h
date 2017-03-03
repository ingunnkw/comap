#ifndef HEALMAPINC
#define HEALMAPINC

#include <map>
#include <sslice.h>
#include <sgrid.h>
#include <serror.h>
#include <stdio.h>

namespace skn {

struct HMapError {};
struct HMapSizeError : public HMapError {};
struct HMapNoBadval  : public HMapError {};

template <typename T> T get_badval() { throw HMapNoBadval(); }
template <> int get_badval<int>() { return -1; }
template <> float get_badval<float>() { return -1.6375e30; }
template <> double get_badval<double>() { return -1.6375e30; }

template <typename T> class HMap;
template <typename T> HMap<T> concatenate(const std::vector<HMap<T> > & maps, int direction = 0);

template <typename T>
class HMap
{
public:
	HMap():nside(1), ordering(0), fullsky(1), data(1,1,12) {}
	HMap(int nside_, int ordering_, int ncomp = 1, int nmap = 1):
		nside(nside_), ordering(ordering_), fullsky(1),
		data(nmap,ncomp,12*nside_*nside_) {}
	HMap(int nside_, int ordering_, const Grid<T,3> & data_):
		nside(nside_), ordering(ordering_), fullsky(1),
		data(data_) { if(data.size(2) != 12*nside*nside) throw HMapSizeError(); }
	HMap(int nside_, int ordering_, const Grid<int,1> pixels_, int ncomp = 1, int nmap = 1):
		nside(nside_), ordering(ordering_), fullsky(0),
		pixels(pixels_), data(nmap,ncomp,pixels_.size()) {}
	HMap(int nside_, int ordering_, const Grid<int,1> pixels_, const Grid<T,3> & data_):
		nside(nside_), ordering(ordering_), fullsky(0),
		pixels(pixels_), data(data_)
		{ if(data.size(2) != pixels.size()) throw HMapSizeError(); }

	// These operations apply to the actual raw data. For a full map,
	// this is the same as the pixels. For a sparse map, it isn't.
	// You may want to use pix() instead.
	T & operator()(int a, int b, int c) { return data(a,b,c); }
	const T & operator()(int a, int b, int c) const { return data(a,b,c); }
	HMap operator()(const Slice & slice) const
	{ return slice_helper(slice, false); }

	T & operator[](int i) { return data[i]; }
	const T & operator[](int i) const { return data[i]; }

	// The apply to the full pixel space, which may have missing values
	T & pix(int a, int b, int c) {
		if(fullsky) return data(a,b,c);
		if(pix2ind.empty()) build(&pix2ind);
		std::map<int,int>::const_iterator i = pix2ind.find(c);
		return i != pix2ind.end() ? data(a,b,i->second) : (tmp=get_badval<T>());
	}
	const T & pix(int a, int b, int c) const {
		static const T tmp = get_badval<T>();
		if(fullsky) return data(a,b,c);
		if(pix2ind.empty()) build(&pix2ind);
		std::map<int,int>::const_iterator i = pix2ind.find(c);
		return i != pix2ind.end() ? data(a,b,i->second) : tmp;
	}
	HMap pix(const Slice & slice) const
	{ return slice_helper(slice, true); }

	int get_nside() const { return nside; }
	int get_ordering() const { return ordering; }
	const Grid<int,1> get_pixels() const { return pixels; } // returns empty if fullsky
	const Vector<int,3> extent() const { return data.extent(); }
	int size(int dim) const { return data.size(dim); }
	int size() const { return data.size();; }
	bool is_fullsky() const { return fullsky; }
	friend HMap concatenate<>(const std::vector<HMap> &, int);

private:
	void build(const std::map<int,int> * p2i) const
	{
		std::map<int,int> * p2iw = const_cast<std::map<int,int>*>(p2i);
		for(int i = 0; i < pixels.size(); i++) (*p2iw)[pixels[i]] = i;
	}
	HMap slice_helper(const Slice & slice, bool pix_slice) const
	{
		int npix = 12*nside*nside;
		int n    = pix_slice ? npix : data.size(2);

		// Setup slice
		std::vector<int> ext(data.extent().size());
		for(int i = 0; i < ext.size(); i++) ext[i] = data.size(i);
		ext[2] = n;
		Slice islice = apply_slice(slice,ext);
		Vector<int,3> oext(islice[0][1],islice[1][1],islice[2][1]);

		/* Copy over data */
		HMap omap;
		omap.nside = nside;
		omap.ordering = ordering;
		omap.fullsky = oext(2) == npix;
		omap.data.resize(oext);
		for(int m = 0, im = islice[0][0]; m < oext(0); m++, im += islice[0][2])
		for(int c = 0, ic = islice[1][0]; c < oext(1); c++, ic += islice[1][2])
		for(int p = 0, ip = islice[2][0]; p < oext(2); p++, ip += islice[2][2])
			omap.data(m,c,p) = pix_slice ? pix(im,ic,ip) : data(im,ic,ip);
		if(oext(2) < npix)
		{
			omap.pixels.resize(oext(2));
			if(pix_slice || pixels.size() == 0)
				for(int p = 0, ip = islice[2][0]; p < oext(2); p++, ip += islice[2][2])
					omap.pixels(p) = ip;
			else
				for(int p = 0, ip = islice[2][0]; p < oext(2); p++, ip += islice[2][2])
					omap.pixels(p) = pixels(ip);
		}
		return omap;
	}

	Grid<T,3> data;
	Grid<int,1> pixels;
	std::map<int,int> pix2ind;
	T tmp;
	int nside, ordering, fullsky;
};


/* Concatenate maps in the submap direction. All parameters must be
 * consistent */
template <typename T>
HMap<T> concatenate(const std::vector<HMap<T> > & maps, int direction)
{
	if(maps.empty()) return HMap<T>();
	HMap<T> res;
	res.nside    = maps[0].nside;
	res.ordering = maps[0].ordering;
	res.fullsky  = maps[0].fullsky;
	res.pixels   = maps[0].pixels;
	// Check consistency and count number of submaps
	Vector<int,3> ext = maps[0].data.extent();
	for(int i = 1; i < maps.size(); i++)
	{
		if(maps[i].nside != maps[0].nside)
			serror("Map %d and %d have inconsistent nside (%d vs %d)!",
				0,i,maps[0].nside, maps[i].nside);
		if(maps[i].ordering != maps[0].ordering)
			serror("Map %d and %d have inconsistent ordering (%d vs %d)!",
				0,i,maps[0].ordering, maps[i].ordering);
		Vector<int,3> ext2 = maps[i].data.extent();
		for(int j = 0; j < ext.size(); j++)
			if(j != direction && ext(j) != ext2(j))
				serror("Map %d and %d have inconsistent dim-%d: %d vs %d!",
					0,i,ext(j),ext2(j));
		ext(direction) += ext2(direction);
	}
	res.data.resize(ext);
	// And copy over. This could probably be more elegant
	if(direction == 0)
	{
		for(int i = 0, j = 0; i < maps.size(); i++)
		for(int k = 0; k < maps[i].data.size(0); k++, j++)
		for(int c = 0; c < maps[i].data.size(1); c++)
		for(int p = 0; p < maps[i].data.size(2); p++)
			res.data(j,c,p) = maps[i].data(k,c,p);
	}
	else if(direction == 1)
	{
		for(int i = 0, j = 0; i < maps.size(); i++)
		for(int c = 0; c < maps[i].data.size(1); c++, j++)
		for(int k = 0; k < maps[i].data.size(0); k++)
		for(int p = 0; p < maps[i].data.size(2); p++)
			res.data(k,j,p) = maps[i].data(k,c,p);
	}
	else if(direction == 2)
	{
		for(int i = 0, j = 0; i < maps.size(); i++)
		for(int p = 0; p < maps[i].data.size(2); p++,j++)
		for(int k = 0; k < maps[i].data.size(0); k++)
		for(int c = 0; c < maps[i].data.size(1); c++)
			res.data(k,c,j) = maps[i].data(k,c,p);
	}
	return res;
}

}

#endif
