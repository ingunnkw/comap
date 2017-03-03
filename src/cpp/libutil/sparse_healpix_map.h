#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <healpix_base.h>
#include <arr.h>
#include <pset.h>
#include <limits.h>
#include <cmath>

 /***************************************************
 * A slow sparse Healpix std::map (but non-adaptive)*
 * Implements the same interface as Healpix_Map     *
 * in order to be easy to use as a replacement.     *
 ***************************************************/
namespace {
template <class T, class S> T Add(const T & a, const S & b) { return a+b; }
template <class T, class S> T Subtract(const T & a, const S & b) { return a-b; }
template <class T, class S> T Multiply(const T & a, const S & b) { return a*b; }
template <class T, class S> T Divide(const T & a, const S & b) { return a/b; }
template <class T, class S> T Power(const T & a, const S & b) { return pow(a,b); }
template <class T, class S> T First(const T & a, const S & b) { return a; }
template <class T, class S> T Second(const T & a, const S & b) { return b; }
template <class T, class S> T Average(const T & a, const S & b) { return (a+b)/2; }
}

typedef Pset<int> Pixset;

template <class T>
class Sparse_Healpix_Map : public Healpix_Base
{
public:
	Sparse_Healpix_Map():nullval(Healpix_undef) {}
	// If no sparsity is passed, then a full std::map is assumed.
	Sparse_Healpix_Map(int order, Healpix_Ordering_Scheme scheme, const std::vector<T> & dat = std::vector<T>()):Healpix_Base(order, scheme),nullval(Healpix_undef) {
		ind2map.resize(Npix());
		for(int i = 0; i < size(); i++) ind2map[i] = i;
		update_map2ind();
		data.resize(ind2map.size());
		for(int i = 0; i < ind2map.size() && i < dat.size(); i++)
			data[i] = dat[i];
	}
	Sparse_Healpix_Map(int nside, Healpix_Ordering_Scheme scheme, const nside_dummy, const std::vector<T> & dat = std::vector<T>()):Healpix_Base(nside, scheme, SET_NSIDE),nullval(Healpix_undef) {
		ind2map.resize(Npix());
		for(int i = 0; i < size(); i++) ind2map[i] = i;
		update_map2ind();
		data.resize(ind2map.size());
		for(int i = 0; i < ind2map.size() && i < dat.size(); i++)
			data[i] = dat[i];
	}
	Sparse_Healpix_Map(int order, const std::vector<int> & sparsity,  Healpix_Ordering_Scheme scheme, const std::vector<T>&dat=std::vector<T>()):Healpix_Base(order, scheme),ind2map(sparsity),nullval(Healpix_undef) {
		update_map2ind();
		data.resize(ind2map.size());
		for(int i = 0; i < ind2map.size() && i < dat.size(); i++)
			data[i] = dat[i];
	}
	Sparse_Healpix_Map(int nside, const std::vector<int> & sparsity,  Healpix_Ordering_Scheme scheme,const nside_dummy, const std::vector<T> & dat = std::vector<T>()):Healpix_Base(nside, scheme, SET_NSIDE),ind2map(sparsity),nullval(Healpix_undef) {
		update_map2ind();
		data.resize(ind2map.size());
		for(int i = 0; i < ind2map.size() && i < dat.size(); i++)
			data[i] = dat[i];
	}
	Sparse_Healpix_Map(const Sparse_Healpix_Map & m):Healpix_Base(m.Order(), m.Scheme()),ind2map(m.ind2map),data(m.data),nullval(m.nullval)
	{ update_map2ind(); }
	// Conversion between different templates
	template <class S>
	Sparse_Healpix_Map(const Sparse_Healpix_Map<S> & m):Healpix_Base(m.Order(),m.Scheme()),ind2map(m.get_sparsity()),nullval(m.get_nullval()),data(m.size())
	{
		for(int i = 0; i < size(); i++) data[i] = m(i);
		update_map2ind();
	}

	// Sets Nside and ordering, and upgrades/degrades the mapping so that
	// it still makes sense
#if 0
	void SetNside(int nside) { SetNside(nside,Scheme()); }
	void SetNside(int nside, Healpix_Ordering_Scheme scheme)
	{
		Healpix_Base other(nside, scheme, SET_NSIDE);
		pix2xyf to_xyf   =       Scheme() == RING ? &Sparse_Healpix_Map::ring2xyf : &Sparse_Healpix_Map::nest2xyf;
		xyf2pix from_xyf = other.Scheme() == RING ? &Sparse_Healpix_Map::xyf2ring : &Sparse_Healpix_Map::xyf2nest;
		int step = other.Order() - Order();
		std::vector<int> i2m;
		std::vector<T> dat;
		if(step >= 0)
		{
			for(int i = 0; i < ind2map.size(); i++)
			{
				int x1, y1, f; (this->*to_xyf)(ind2map[i], x1, y1, f);
				x1 <<= step; y1 <<= step;
				int x2 = x1 + (1 << step), y2 = y1 + (1 << step);
				for(int x = x1; x < x2; x++)
				for(int y = y1; y < y2; y++)
				{
					int opix = (other.*from_xyf)(x,y,f);
					i2m.push_back(opix);
					dat.push_back(data[i]);
				}
			}
		}
		else
		{
			std::map<int,int> counts;
			std::map<int,T> vals;
			for(int i = 0; i < ind2map.size(); i++)
			{
				int x, y, f; (this->*to_xyf)(ind2map[i], x, y, f);
				int opix = (other.*from_xyf)(x>>-step,y>>-step,f);
				counts[opix]++;
				vals[opix] += data[i];
			}
			std::map<int,int>::const_iterator i;
			typename std::map<int,T>::const_iterator j;
			for(i = counts.begin(), j = vals.begin(); i != counts.end() && j != vals.end(); i++, j++)
			{
				dat.push_back(j->second/i->second);
				i2m.push_back(i->first);
			}
		}
		Healpix_Base::SetNside(nside, scheme);
		data = dat;
		ind2map = i2m;
		update_map2ind();
	}
	void SetNside(int nside, const std::vector<int> & sparsity) { SetNside(nside,sparsity, Scheme()); }
	void SetNside(int nside, const std::vector<int> & sparsity, Healpix_Ordering_Scheme scheme)
	{
		Sparse_Healpix_Map<T> tmp(*this);
		Healpix_Base::SetNside(nside, scheme);
		ind2map = sparsity;
		data.resize(sparsity.size());
		update_map2ind();
		import(tmp);
	}
#endif
	void fill(const T & val) { for(int i = 0; i < data.size(); i++) data[i]=val; }
#if 0
	void import(const Sparse_Healpix_Map & other, bool pessimistic = false)
	{
		int step = other.Order() - Order();
		if(step >= 0)
		{
			pix2xyf to_xyf   =       Scheme() == RING ? &Sparse_Healpix_Map::ring2xyf : &Sparse_Healpix_Map::nest2xyf;
			xyf2pix from_xyf = other.Scheme() == RING ? &Sparse_Healpix_Map::xyf2ring : &Sparse_Healpix_Map::xyf2nest;
			for(int i = 0; i < ind2map.size(); i++)
			{
				int pix = ind2map[i];
				int x1, y1, f; (this->*to_xyf)(pix, x1, y1, f);
				x1 <<= step; y1 <<= step;
				int x2 = x1 + 1 << step, y2 = y1 + 1 << step;
				int n = 0;
				T val = 0;
				bool ok = true;
				for(int x = x1; x < x2 && ok; x++)
				for(int y = y1; y < y2 && ok; y++)
				{
					int opix = (other.*from_xyf)(x,y,f);
					std::map<int,int>::const_iterator it = other.map2ind.find(opix);
					if(it != other.map2ind.end()) {
						val += other.data[it->second];
						n++;
					}
					else if(pessimistic) ok = false;
				}
				ok &= n != 0;
				data[i] = ok ? val/n : nullval;
			}
		}
		else
		{
			pix2xyf to_xyf   = other.Scheme() == RING ? &Sparse_Healpix_Map::ring2xyf : &Sparse_Healpix_Map::nest2xyf;
			xyf2pix from_xyf =       Scheme() == RING ? &Sparse_Healpix_Map::xyf2ring : &Sparse_Healpix_Map::xyf2nest;
			// basically as above, but the other way around
			for(int i = 0; i < other.ind2map.size(); i++)
			{
				int opix = other.ind2map[i];
				int x1, y1, f; (other.*to_xyf)(opix, x1, y1, f);
				x1 >>= step; y1 >>= step;
				int x2 = x1 + 1 >> step, y2 = y1 + 1 >> step;
				T val = other.data[i];
				for(int x = x1; x < x2; x++)
				for(int y = y1; y < y2; y++)
				{
					int pix = (this->*from_xyf)(x,y,f);
					std::map<int,int>::const_iterator it = map2ind.find(pix);
					if(it != map2ind.end()) data[it->second] = val;
				}
			}
		}
	}
	// These just call import. They are provided for compatibility
	void import_upgrade(const Sparse_Healpix_Map & other) { import(other); }
	void import_nograde(const Sparse_Healpix_Map & other) { import(other); }
	void import_degrade(const Sparse_Healpix_Map & other, bool pessimistic = false) { import(other, pessimistic); }
#endif
	T & operator[](int i) {
		std::map<int,int>::const_iterator j = map2ind.find(i);
		return j != map2ind.end() ? data[j->second] : scratch = nullval;
	}
	const T & operator[](int i) const {
		std::map<int,int>::const_iterator j = map2ind.find(i);
		return j != map2ind.end() ? data[j->second] : nullval;
	}
	T & operator()(int i) { return data[i]; }
	const T & operator()(int i) const { return data[i]; }
	int size() const { return ind2map.size(); }
	int index_to_pixel(int i) const { return ind2map[i]; }
	int pixel_to_index(int i) const {
		std::map<int,int>::const_iterator j = map2ind.find(i);
		return j != map2ind.end() ? j->second : -1;
	}
	bool defined(int i) const { return map2ind.find(i) != map2ind.end(); }
	void swap_scheme() {
		if(Scheme() == RING)
		{
			for(int i = 0; i < ind2map.size(); i++)
				ind2map[i] = ring2nest(ind2map[i]);
			scheme_ = NEST;
		}
		else
		{
			for(int i = 0; i < ind2map.size(); i++)
				ind2map[i] = nest2ring(ind2map[i]);
			scheme_ = RING;
		}
		update_map2ind();
	}
	T interpolation(const fix_arr<int,4> & pix, const fix_arr<double,4> & wgt) const
	{
		int n = 0; double wsum = 0; T res = 0;
		for(int i = 0; i < 4; i++)
		{
			std::map<int,int>::const_iterator it = map2ind.find(pix[i]);
			if(it != map2ind.end()) {
				n++; wsum += wgt[i]; res += data[it->second]*wgt[i];
			}
			return n ? res/wsum : nullval;
		}
		return T();
	}
	T interpolated_value(const pointing & ptg) const
	{
		fix_arr<int,4> pix; fix_arr<double,4> wgt;
		get_interpol(ptg,pix,wgt);
		return interpolation(pix,wgt);
	}
	void swap(Sparse_Healpix_Map<T> & other) {
		data.swap(other.data);
		ind2map.swap(other.ind2map);
		map2ind.swap(other.map2ind);
		Healpix_Base::swap(other);
	}
	void minmax(T&Min,T&Max) const {
		Min = data[0],Max=data[0];
		for(int i = 1; i <size(); i++)
		{
			if(data[i] < Min) Min = data[i];
			if(data[i] > Max) Max = data[i];
		}
	}
	T sum() const {
		if(size() == 0) return nullval;
		T res();
		for(int i = 0; i < size(); i++) res += data[i];
		return res;
	}
	T average() const { return sum()/size(); }
	T rms() const {
		if(size() == 0) return nullval;
		T mean = average(), res();
		for(int i = 0; i < size(); i++) { T t = data[i]-mean; res += t*t; }
		return sqrt(res/size());
	}
	T absmax() const {
		if(size() == 0) return nullval;
		T m = abs(data[0]);
		for(int i = 1; i < size(); i++) { T a = abs(data[i]); if(a > m) m = a; }
		return m;
	}
	void add(const T & val) { for(int i = 0; i < size(); i++) data[i]+=val; }
	void set_nullval(const T & t) { nullval = t; }
	const T & get_nullval() const { return nullval; }
	void set_sparsity(const std::vector<int> & v) {
		std::vector<T> dat(v.size());
		for(int i = 0; i < v.size(); i++)
			dat[i] = (*this)[v[i]];
		ind2map = v;
		data = dat;
		update_map2ind();
	}
	const std::vector<int> & get_sparsity() const { return ind2map; }
	void prune() { prune(nullval*1.00001,nullval/1.00001); }
	void prune(const T & v) { prune(v,v); }
	void prune(const T & v1, const T & v2)
	{
		std::vector<double> v;
		std::vector<int> ind;
		bool arg_nan = v1 != v1 || v2 != v2;
		for(int i = 0; i < size(); i++)
		{
			bool is_nan = data[i] != data[i];
			if(arg_nan && is_nan) continue;
			if(data[i] >= v1 && data[i] <= v2) continue;
			v.push_back(data[i]);
			ind.push_back(ind2map[i]);
		}
		ind2map = ind;
		data = v;
		update_map2ind();
	}

	// The above is mostly what is in Healpix_Map. Now for some stuff
	// that should have been there:
	
	// Operators. Maps should be the same ordering before you try any of these.
	// Not necessarily the same sparsity, though. The result is the intersection,
	// Except for the | operator, which is the union
	Sparse_Healpix_Map & operator=(const Sparse_Healpix_Map & m)
	{
		if(this == &m) return *this;
		Healpix_Base::SetNside(m.Nside(),m.Scheme());
		data = m.data;
		ind2map = m.ind2map;
		update_map2ind();
		return *this;
	}
	template <class S>
	Sparse_Healpix_Map<T> & operator+=(const Sparse_Healpix_Map<S> & m) { return *this = operate(m,Add<T,S>); }
	template <class S>
	Sparse_Healpix_Map<T> & operator-=(const Sparse_Healpix_Map<S> & m) { return *this = operate(m,Subtract<T,S>); }
	template <class S>
	Sparse_Healpix_Map<T> & operator*=(const Sparse_Healpix_Map<S> & m) { return *this = operate(m,Multiply<T,S>); }
	template <class S>
	Sparse_Healpix_Map<T> & operator/=(const Sparse_Healpix_Map<S> & m) { return *this = operate(m,Divide<T,S>); }
	template <class S>
	Sparse_Healpix_Map<T> operator&=(const Sparse_Healpix_Map<S> & m) const { return *this = operate(m,First<T,S>); }
	template <class S>
	Sparse_Healpix_Map<T> operator+(const Sparse_Healpix_Map<S> & m) const { return operate(m,Add<T,S>); }
	template <class S>
	Sparse_Healpix_Map<T> operator-(const Sparse_Healpix_Map<S> & m) const { return operate(m,Subtract<T,S>); }
	template <class S>
	Sparse_Healpix_Map<T> operator*(const Sparse_Healpix_Map<S> & m) const { return operate(m,Multiply<T,S>); }
	template <class S>
	Sparse_Healpix_Map<T> operator/(const Sparse_Healpix_Map<S> & m) const { return operate(m,Divide<T,S>); }
	template <class S>
	Sparse_Healpix_Map<T> pow(const Sparse_Healpix_Map<S> & m) const { Sparse_Healpix_Map res(*this); res.operate(m,Power<T,S>); return res; }
	template <class S>
	Sparse_Healpix_Map<T> operator&(const Sparse_Healpix_Map<S> & m) const { return operate(m,First<T,S>); }
	template <class S>
	Sparse_Healpix_Map<T> operator|(const Sparse_Healpix_Map<S> & m) const
	{
		std::map<int,int> counts;
		std::map<int,T> vals;
		for(int i = 0; i < size(); i++)
		{
			counts[ind2map[i]]++;
			vals[ind2map[i]] += data[i];
		}
		for(int i = 0; i < m.size(); i++)
		{
			counts[m.ind2map[i]]++;
			vals[m.ind2map[i]] += m.data[i];
		}
		std::vector<int> pix;
		std::vector<T> dat;
		std::map<int,int>::const_iterator i;
		typename std::map<int,T>::const_iterator j;
		for(i = counts.begin(), j = vals.begin(); i != counts.end() && j != vals.end(); i++, j++)
		{ dat.push_back(j->second/i->second); pix.push_back(i->first); }
		return Sparse_Healpix_Map(Order(),pix,Scheme(),dat);
	}

	template <class S>
	Sparse_Healpix_Map<T> & operator=(const S & m) { operate(m,Second<T,S>); return*this;}
	template <class S>
	Sparse_Healpix_Map<T> & operator+=(const S & m) { operate(m,Add<T,S>); return*this;}
	template <class S>
	Sparse_Healpix_Map<T> & operator-=(const S & m) { operate(m,Subtract<T,S>);return*this; }
	template <class S>
	Sparse_Healpix_Map<T> & operator*=(const S & m) { operate(m,Multiply<T,S>);return*this; }
	template <class S>
	Sparse_Healpix_Map<T> & operator/=(const S & m) { operate(m,Divide<T,S>);return*this; }
	template <class S>
	Sparse_Healpix_Map<T> operator+(const S & m) const { Sparse_Healpix_Map res(*this); res.operate(m,Add<T,S>); return res;}
	template <class S>
	Sparse_Healpix_Map<T> operator-(const S & m) const { Sparse_Healpix_Map res(*this); res.operate(m,Subtract<T,S>);return res; }
	template <class S>
	Sparse_Healpix_Map<T> operator*(const S & m) const { Sparse_Healpix_Map res(*this); res.operate(m,Multiply<T,S>);return res; }
	template <class S>
	Sparse_Healpix_Map<T> operator/(const S & m) const { Sparse_Healpix_Map res(*this); res.operate(m,Divide<T,S>);return res; }
	template <class S>
	Sparse_Healpix_Map<T> pow(const S & m) const { Sparse_Healpix_Map res(*this); res.operate(m,Power<T,S>); return res; }
	Sparse_Healpix_Map<T> operator-() const { Sparse_Healpix_Map res(*this); for(int i = 0; i < res.size(); i++) res.data[i] = -res.data[i]; return res; }

	std::vector<Pixset> neighbor_number() const
	{
		std::vector<Pixset> res(9);
		for(int i = 0; i < ind2map.size(); i++)
		{
			fix_arr<int,8> neighs;
			neighbors(ind2map[i], neighs);
			int n = 0;
			for(int j = 0; j < neighs.size(); j++)
				n += neighs[j] != -1 && pixel_to_index(neighs[j]) != -1;
			res[n] += ind2map[i];
		}
		return res;
	}

	std::vector<Pixset> neighbors_missing() const
	{
		std::vector<Pixset> res(9);
		for(int i = 0; i < ind2map.size(); i++)
		{
			fix_arr<int,8> neighs;
			neighbors(ind2map[i], neighs);
			int n = 0;
			for(int j = 0; j < neighs.size(); j++)
				n += neighs[j] != -1 && pixel_to_index(neighs[j]) == -1;
			res[n] += ind2map[i];
		}
		return res;
	}

	// A pretty heavy, but sometimes needed operation
	std::map<int,Pixset> layers(int max = INT_MAX, int min = 0) const
	{
		std::map<int,Pixset> res;
		Pixset & border = res[0];
		// First find boundary pixels
		std::vector<Pixset> nneigh = neighbors_missing();
		for(int i = 1; i < nneigh.size(); i++)
			border += nneigh[i];
		Pixset inside = ind2map;
		// Now iterate inwards
		for(int i = 1; i <= max && !res[i-1].empty(); i++)
		{
			for(Pixset::const_iterator j = res[i-1].begin(); j != res[i-1].end(); j++)
			{
				fix_arr<int,8> neigh; neighbors(*j,neigh);
				for(int k = 0; k < neigh.size(); k++)
					if(!(neigh[k] < 0 || i == 1 && !(neigh[k] < inside) || i > 1 && neigh[k] < res[i-2] || neigh[k] < res[i-1]))
						res[i] += neigh[k];
			}
		}
		// And outwards
		for(int i = -1; i >= min && !res[i+1].empty(); i--)
		{
			for(Pixset::const_iterator j = res[i+1].begin(); j != res[i+1].end(); j++)
			{
				fix_arr<int,8> neigh; neighbors(*j,neigh);
				for(int k = 0; k < neigh.size(); k++)
					if(!(neigh[k] < 0 || i == -1 && neigh[k] < inside || i < -1 && neigh[k] < res[i+2] || neigh[k] < res[i+1]))
						res[i] += neigh[k];
			}
		}
		return res;
	}
	double pixsize() const { return sqrt(4*M_PI/Npix()); }
protected:
	template <class S>
	Sparse_Healpix_Map<T> operate(const Sparse_Healpix_Map<S> & m, T (*operation)(const T&,const S&)) const
	{
		std::vector<int> pix;
		std::vector<T> dat;
		for(int i = 0; i < size(); i++)
		{
			std::map<int,int>::const_iterator j = m.map2ind.find(ind2map[i]);
			if(j == m.map2ind.end()) continue;
			pix.push_back(j->first);
			dat.push_back(operation(data[i],m.data[j->second]));
		}
		return Sparse_Healpix_Map<T>(Order(),pix,Scheme(),dat);
	}
	// This changes the current std::map
	template <class S>
	void operate(const S & val, T (*operation)(const T &, const S&))
	{ for(int i = 0; i < size(); i++) data[i] = operation(data[i],val); }
	void update_map2ind() {
		map2ind.clear();
		for(int i = 0; i < ind2map.size(); i++)
			map2ind[ind2map[i]] = i;
	}
	int nside2order(int nside) { int i = -1; for(;nside;nside >>= 1, i++); return i; }
	T nullval, scratch;
	std::vector<T> data;
	std::vector<int> ind2map;
	std::map<int,int> map2ind;
};

template <class T> Sparse_Healpix_Map<T> pow(const Sparse_Healpix_Map<T> & a, const Sparse_Healpix_Map<T> & b)
{ return a.pow(b); }
template <class T> Sparse_Healpix_Map<T> pow(const Sparse_Healpix_Map<T> & a, const T & b)
{ return a.pow(b); }

template <class T>
Sparse_Healpix_Map<T> operator+(const T & m, const Sparse_Healpix_Map<T> & n) { return n+m; }
template <class T>
Sparse_Healpix_Map<T> operator-(const T & m, const Sparse_Healpix_Map<T> & n)
{
	Sparse_Healpix_Map<T> res(n);
	for(int i = 0; i < res.size(); i++) res(i) = m-res(i);
	return res;
}
template <class T>
Sparse_Healpix_Map<T> operator*(const T & m, const Sparse_Healpix_Map<T> & n) { return n*m; }
template <class T>
Sparse_Healpix_Map<T> operator/(const T & m, const Sparse_Healpix_Map<T> & n)
{
	Sparse_Healpix_Map<T> res(n);
	for(int i = 0; i < res.size(); i++) res(i) = m/res(i);
	return res;
}
template <class T>
Sparse_Healpix_Map<T> pow(const T & m, const Sparse_Healpix_Map<T> & n)
{
	Sparse_Healpix_Map<T> res(n);
	for(int i = 0; i < res.size(); i++) res(i) = std::pow(m, n(i));
	return res;
}

// Put all maps on common sparsity by expanding the sparsity of each map
template<class T>
void unify(std::vector<Sparse_Healpix_Map<T> > & maps)
{
	std::set<int> pixels;
	for(int i = 0; i < maps.size(); i++)
		for(int j = 0; j < maps[i].size(); j++)
			pixels.insert(maps[i].index_to_pixel(j));
	std::vector<int> pix;
	for(std::set<int>::const_iterator i = pixels.begin(); i != pixels.end();i++)
		pix.push_back(*i);
	for(int i = 0; i < maps.size(); i++)
		maps[i].set_sparsity(pix);
}

template<class T>
void prune(std::vector<Sparse_Healpix_Map<T> > & maps) { for(int i = 0; i < maps.size(); i++) maps[i].prune(); }
template<class T>
void prune(std::vector<Sparse_Healpix_Map<T> > & maps, const T & v) { for(int i = 0; i < maps.size(); i++) maps[i].prune(v); }
template<class T>
void prune(std::vector<Sparse_Healpix_Map<T> > & maps, const T & v1, const T & v2) { for(int i = 0; i < maps.size(); i++) maps[i].prune(v1,v2); }
