#ifndef FITSTYPESH
#define FITSTYPESH

#include <serror.h>
#include <fitsio.h>

namespace skn {

template<class T> int fitstype() { serror("Fitstype undefined for this type!"); return 0; }
template<> int fitstype<unsigned char>() { return PLANCK_UINT8; }
template<> int fitstype<bool>()          { return PLANCK_BOOL; }
template<> int fitstype<signed short>()  { return PLANCK_INT16; }
template<> int fitstype<signed int>()    { return PLANCK_INT32; }
template<> int fitstype<float>()         { return PLANCK_FLOAT32; }
template<> int fitstype<double>()        { return PLANCK_FLOAT64; }
template<> int fitstype<long long>()     { return PLANCK_INT64; }

int dtype2fits(const std::string & dtype)
{
	if(dtype == "double") return fitstype<double>();
	if(dtype == "float")  return fitstype<float>();
	if(dtype == "int")    return fitstype<int>();
	if(dtype == "byte")   return fitstype<unsigned char>();
	if(dtype == "long")   return fitstype<long long>();
	if(dtype == "bool")   return fitstype<bool>();
	serror("Unknown data type %s", dtype.c_str());
	return fitstype<double>();
}

}

#endif
