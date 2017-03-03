#ifndef AGRIDA
#define AGRIDA

#include <svector.h>
#include <cmath>
#include <algorithm>

namespace skn {

/* The previous Grid class became too hacky. Again, I choose macros as
 * the lesser evil. Where Vec is specialized for short mathematical
 * vectors, this class is slower but more general. */

#define GRIDCOMMON(N) \
	typedef Vector<int,N> Index; \
	Grid(const Index & dim = Index(1)) { setup(dim); } \
	Grid(const Index & dim, const T & t) { setup(dim); set(t); } \
	Grid(const Index & dim, const T * t) { setup(dim); set(t); } \
	Grid(const Grid<T,N> & rhs):l(rhs.l),n(rhs.n),c(rhs.c) { \
		d = new T[l]; \
		for(int i = 0; i < l; i++) d[i] = rhs.d[i]; \
	} \
	~Grid() { delete [] d; } \
	T & operator()(const Index & ind) { return d[lookup(ind)]; } \
	const T & operator()(const Index & ind) const{ return d[lookup(ind)]; } \
	T & operator [](int i) { return d[i]; } \
	const T & operator [](int i) const { return d[i]; } \
	int size() const { return l; } \
	int size(int i) const { return n[i]; } \
	const Index & extent() const { return n; } \
	void resize(const Index & dim) { if(dim == n) return; delete [] d; setup(dim); } \
	template <class S> \
	void set(const S & val) { for(int i = 0; i < size(); i++) d[i] = val; } \
	template <class S> \
	void set(const S * val) { for(int i = 0; i < size(); i++) d[i] = val[i]; } \
	template <class S> \
	Grid<T,N> & operator=(const S & t) { set(t); return *this; } \
	Grid<T,N> & operator=(const Grid<T,N> & rhs) { \
		if(&rhs != this) \
		{ \
			l = rhs.l; n = rhs.n; c = rhs.c; \
			delete [] d; \
			d = new T[l]; \
			for(int i = 0; i < size(); i++) d[i] = rhs.d[i]; \
		} \
		return *this; \
	} \
	void swap(Grid & rhs) { std::swap(l, rhs.l); std::swap(n, rhs.n); std::swap(c, rhs.c); std::swap(d, rhs.d); } \
	int index2offset(const Index & ind) const { return lookup(ind); } \
	Index offset2index(int i) const { \
		Index res; \
		for(int j = 0; j < N-1; j++) \
		{ \
			res[j] = i/c[j]; \
			i -= res[j]*c[j]; \
		} \
		res[N-1] = i; \
		return res; \
	} \
private: \
	void setup(const Index & dim) \
	{ \
		n = dim; \
		if(N>0) { \
			c[N-1] = 1; \
			for(int i = N-2; i >= 0; i--) c[i] = c[i+1]*n[i+1]; \
			l = c[0]*n[0]; \
		} else { l = 1; } \
		d = new T[l]; \
	} \
	int lookup(const Index & dim) const \
	{ \
		int j = 0; \
		for(int i = 0; i < N; i++) j += dim[i]*c[i]; \
		return j; \
	} \
	int l; \
	Index n, c; \
	T * d;

template <class T, int N> class Grid
{
public:
	GRIDCOMMON(N);
};

template <class T> class Grid<T,1>
{
public:
	explicit Grid(int n1, const T & t = T()) { setup(Index(n1)); set(t); }
	Grid(int n1, const T * t) { setup(Index(n1)); set(t); }
	T & operator()(int a) { return d[a]; }
	const T & operator()(int a) const { return d[a]; }
	void resize(int a) { setup(Index(a)); }
	GRIDCOMMON(1);
};

template <class T> class Grid<T,2>
{
public:
	Grid(int n1, int n2, const T & t = T()) { setup(Index(n1,n2)); set(t); }
	Grid(int n1, int n2, const T * t) { setup(Index(n1,n2)); set(t); }
	T & operator()(int a, int b) { return d[a*c[0]+b]; }
	const T & operator()(int a, int b) const { return d[a*c[0]+b]; }
	void resize(int a, int b) { setup(Index(a,b)); }
	GRIDCOMMON(2);
};

template <class T> class Grid<T,3>
{
public:
	Grid(int n1, int n2, int n3, const T & t = T()) { setup(Index(n1,n2,n3)); set(t); }
	Grid(int n1, int n2, int n3, const T * t) { setup(Index(n1,n2,n3)); set(t); }
	T & operator()(int n1, int n2, int n3) { return d[n1*c[0]+n2*c[1]+n3]; }
	const T & operator()(int n1, int n2, int n3) const { return d[n1*c[0]+n2*c[1]+n3]; }
	void resize(int n1, int n2, int n3) { setup(Index(n1,n2,n3)); }
	GRIDCOMMON(3);
};

template <class T> class Grid<T,4>
{
public:
	Grid(int n1, int n2, int n3, int n4, const T & t = T()) { setup(Index(n1,n2,n3,n4)); set(t); }
	Grid(int n1, int n2, int n3, int n4, const T * t) { setup(Index(n1,n2,n3,n4)); set(t); }
	T & operator()(int n1, int n2, int n3, int n4) { return d[lookup(Index(n1,n2,n3,n4))]; }
	const T & operator()(int n1, int n2, int n3, int n4) const { return d[lookup(Index(n1,n2,n3,n4))]; }
	void resize(int n1, int n2, int n3, int n4) { setup(Index(n1,n2,n3,n4)); }
	GRIDCOMMON(4);
};

template <class T> class Grid<T,5>
{
public:
	Grid(int n1, int n2, int n3, int n4, int n5, const T & t = T()) { setup(Index(n1,n2,n3,n4,n5)); set(t); }
	Grid(int n1, int n2, int n3, int n4, int n5, const T * t) { setup(Index(n1,n2,n3,n4,n5)); set(t); }
	T & operator()(int n1, int n2, int n3, int n4, int n5) { return d[lookup(Index(n1,n2,n3,n4,n5))]; }
	const T & operator()(int n1, int n2, int n3, int n4, int n5) const { return d[lookup(Index(n1,n2,n3,n4,n5))]; }
	void resize(int n1, int n2, int n3, int n4, int n5) { setup(Index(n1,n2,n3,n4,n5)); }
	GRIDCOMMON(5);
};

#undef GRIDCOMMON

// The usual operators
#define OPS(OP) \
template <class T, int N> Grid<T,N> & operator OP##=(Grid<T,N> & g, const Grid<T,N> & v) { for(int i = 0; i < g.size(); i++) g[i] OP##= v[i]; return g; } \
template <class T, class S, int N> Grid<T,N> & operator OP##=(Grid<T,N> & g, const S & v) { for(int i = 0; i < g.size(); i++) g[i] OP##= v; return g; } \
template <class T, class S, int N> Grid<T,N> operator OP(const Grid<T,N> & g, const S & v) { Grid<T,N> r(g); return r OP##= v; } \
template <class T, int N> Grid<T,N> operator OP(const Grid<T,N> & g, const Grid<T,N> & v) { Grid<T,N> r(g); return r OP##= v; } \
template <class T, class S, int N> Grid<T,N> operator OP(const S & g, const Grid<T,N> & v) { Grid<T,N> r(v); for(int i = 0; i < r.size(); i++) r[i] = v OP g[i]; return r; }
OPS(+); OPS(-); OPS(*); OPS(/); OPS(|); OPS(&); OPS(^); OPS(<<); OPS(>>);
#undef OPS

#define UNOP(OP) \
template <class T, int N> Grid<T,N> operator OP(const Grid<T,N> & g) { Grid<T,N> r(g); for(int i = 0; i < r.size(); i++) r[i] = OP g[i]; return r; }
UNOP(-); UNOP(~); UNOP(!);
#undef UNOP

// Nifty operations
template <class T, int N> int size(const Grid<T,N> & g) { return g.size(); }
template <class T, int N> int size(const Grid<T,N> & g, int i) { return g.size(i); }
template <class T, int N> T sum(const Grid<T,N> & g) { T res = T(); for(int i = 0; i < g.size(); i++) res += g[i]; return res; }
template <class T, int N> T mean(const Grid<T,N> & g) { return sum(g)/size(g); }
template <class T, int N> int count(const Grid<T,N> & g) { int res = int(); for(int i = 0; i < g.size(); i++) res += g[i]; return res; }
template <class T, int N> const typename Grid<T,N>::Index & extent(const Grid<T,N> & g) { return g.extent(); }
template <class T, int N> void resize(const Grid<T,N> & g, const typename Grid<T,N>::Index & ind) { return g.resize(ind); }
template <class T, int N> T max(const Grid<T,N> & g) { T res = g[0]; for(int i = 1; i < g.size(); i++) res = res < g[i] ? g[i] : res; return res; }
template <class T, int N> T min(const Grid<T,N> & g) { T res = g[0]; for(int i = 1; i < g.size(); i++) res = res < g[i] ? res : g[i]; return res; }

#define MFUN(FUN) \
template <class T, int N> void FUN(const Grid<T,N> & g, Grid<T,N> & res) { for(int i = 0; i < g.size(); i++) res[i] = FUN(g[i]); } \
template <class T, int N> Grid<T,N> FUN(const Grid<T,N> & g) { Grid<T,N> res(g.extent()); FUN(g, res); return res; }
MFUN(sin); MFUN(cos); MFUN(log); MFUN(exp); MFUN(tan); MFUN(asin); MFUN(acos); MFUN(atan); MFUN(sqrt); MFUN(floor); MFUN(ceil); MFUN(sinh); MFUN(cosh); MFUN(tanh); MFUN(abs);
#undef MFUN

}

namespace std
{
	template <class T, int N>
	void swap(skn::Grid<T,N> & a, skn::Grid<T,N> & b) { a.swap(b); }
}

#endif
