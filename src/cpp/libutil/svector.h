#ifndef AVECT
#define AVECT

#include <cmath>

namespace skn {

// Class definitions. A macro seemed to be the lesser evil here.
// Inheritance looks like it might work, but doesn't work out.

#define VECCOMMON(Q) \
	explicit Vector(const T & t = T()) { fill(p,t,Q); } \
	template <class S> \
	Vector(const Vector<S,Q> & v) { copy(p,&v(0),Q); } \
	template <class S> \
	Vector<T,Q> & operator=(const S & t) { fill(p,t,Q); return *this; } \
	Vector<T,Q> & operator=(const Vector<T,Q>&v) { copy(p,v.p,Q); return *this; } \
	T & operator()(int i) { return p[i]; } \
	const T & operator()(int i) const { return p[i]; } \
	T & operator[](int i) { return p[i]; } \
	const T & operator[](int i) const { return p[i]; } \
	int size() const { return Q; } \
private: \
	template <class S> static void fill(T * a, const S & b, int n) \
	{ for(int i = 0; i < Q; i++) a[i] = b; } \
	template <class S> static void copy(T * a, const S * b, int n) \
	{ for(int i = 0; i < Q; i++) a[i] = b[i]; } \
	T p[Q];

template <class T, int N> class Vector
{
public:
	VECCOMMON(N);
};

template <class T> class Vector<T,2>
{
public:
	Vector(const T & a, const T & b) { p[0] = a; p[1] = b; }
	VECCOMMON(2);
};

template <class T> class Vector<T,3>
{
public:
	Vector(const T & a, const T & b, const T & c) { p[0] = a; p[1] = b; p[2] = c;}
	VECCOMMON(3);
};

template <class T> class Vector<T,4>
{
public:
	Vector(const T & a, const T & b, const T & c, const T & d) { p[0] = a; p[1] = b; p[2] = c; p[3] = d; }
	VECCOMMON(4);
};

template <class T> class Vector<T,5>
{
public:
	Vector(const T & a, const T & b, const T & c, const T & d, const T & e) { p[0] = a; p[1] = b; p[2] = c; p[3] = d; p[4] = e; }
	VECCOMMON(5);
};

template <class T> class Vector<T,6>
{
public:
	Vector(const T & a, const T & b, const T & c, const T & d, const T & e, const T & f) { p[0] = a; p[1] = b; p[2] = c; p[3] = d; p[4] = e; p[5] = f; }
	VECCOMMON(6);
};

#undef VECCOMMON

// Similar operators. I am starting to see the point of temporary macros.
#define OPS(OP) \
template <class T, int N> Vector<T,N> & operator OP##=(Vector<T,N> & a, const Vector<T,N> & b) { for(int i = 0; i < N; i++) a[i] OP##= b[i]; return a;} \
template <class T, class S, int N> Vector<T,N> & operator OP##=(Vector<T,N> & a, const S & b) { for(int i = 0; i < N; i++) a[i] OP##= b; return a;} \
template <class T, int N> Vector<T,N> operator OP(const Vector<T,N> & a, const Vector<T,N> & b) { Vector<T,N> r(a); return r OP##= b;} \
template <class T, class S, int N> Vector<T,N> operator OP(const Vector<T,N> & a, const S & b) { Vector<T,N> r(a); return r OP##= b;} \
template <class T, class S, int N> Vector<T,N> operator OP(const T & a, const Vector<T,N> & v) { Vector<T,N> r; for(int i = 0; i < N; i++) r[i] = a OP v[i]; return r; }
OPS(+); OPS(-); OPS(*); OPS(/); OPS(|); OPS(&); OPS(^); OPS(%); OPS(<<); OPS(>>);
#undef OPS

#define UNOP(OP) template <class T, int N> Vector<T,N> operator OP(const Vector<T,N> & a) { Vector<T,N> r; for(int i = 0; i < N; i++) r[i] = OP a[i]; return r; }
UNOP(-); UNOP(~);
#undef UNOP

// Comparisons
template <class T, int N>
int compare(const Vector<T,N> & a, const Vector<T,N> & b)
{
	for(int i = 0; i < N; i++)
		if(a[i] > b[i]) return 1;
		else if(a[i] < b[i]) return -1;
	return 0;
}
#define COMPS(OP) template <class T, int N> bool operator OP(const Vector<T,N> & a, const Vector<T,N> & b) { return compare(a,b) OP 0; }
COMPS(==); COMPS(!=); COMPS(>); COMPS(<);
#undef COMPS

// Useful functions
template <class T, int N> T sum(const Vector<T,N> & a) { T res = T(); for(int i = 0; i < N; i++) res += a[i]; return res; }
template <class T, int N> T dot(const Vector<T,N> & a, const Vector<T,N> & b) { T res = T(); for(int i = 0; i < N; i++) res += a[i]*b[i]; return res; }
template <class T, int N> T sqr(const Vector<T,N> & a) { return dot(a,a); }
template <class T, int N> T abs(const Vector<T,N> & a) { return std::sqrt(dot(a,a)); }
template <class T, int N> T len(const Vector<T,N> & a) { return abs(a); }
template <class T, int N> int size(const Vector<T,N> & a) { return N; }
template <class T, int N> T angle(const Vector<T,N> & a, const Vector<T,N> & b) { return std::acos(dot(a,b)/sqrt(sqr(a)*sqr(b))); }
template <class T, int N> Vector<T,N> & rescale(Vector<T,N> & a, const T & l) { return a *= l/len(a); }
template <class T, int N> Vector<T,N> rescaled(const Vector<T,N> & a, const T & l) { Vector<T,N> r(a); return rescale(r, l); }
template <class T, int N> T overlap(const Vector<T,N> & a, const Vector<T,N> & b) { return dot(a,b)/len(a); }
template <class T, int N> Vector<T,N> proj(const Vector<T,N> & a, const Vector<T,N> & b) { return rescaled(a, overlap(a,b)); }
template <class T, int N> int argmax(const Vector<T,N> & a) {
	int i,j; for(i = 1, j = 0; i < N; i++) if(a[i] > a[j]) j = i;
	return j;
}
template <class T, int N> int argmin(const Vector<T,N> & a) {
	int i,j; for(i = 1, j = 0; i < N; i++) if(a[i] < a[j]) j = i;
	return j;
}
template <class T, int N> const T & max(const Vector<T,N> & a) { return a[argmax(a)]; }
template <class T, int N> const T & min(const Vector<T,N> & a) { return a[argmin(a)]; }

#define MFUN(FUN) \
template <class T, int N> Vector<T,N> FUN(const Vector<T,N> & a) { Vector<T,N> res; for(int i = 0; i < N; i++) res[i] = FUN(a[i]); }
MFUN(sin); MFUN(cos); MFUN(log); MFUN(exp); MFUN(tan); MFUN(asin); MFUN(acos); MFUN(atan); MFUN(sqrt); MFUN(floor); MFUN(ceil); MFUN(sinh); MFUN(cosh); MFUN(tanh);
#undef MFUN

// The specialized functions:
template <class T> Vector<T,2> hat(const Vector<T,2> & a) { return Vector<T,2>(-a[1],a[0]); }
template <class T> Vector<T,3> cross(const Vector<T,3> & a, const Vector<T,3> & b) { return Vector<T,3>(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]); }

}

#endif
