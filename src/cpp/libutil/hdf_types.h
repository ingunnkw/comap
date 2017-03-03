#ifndef HDFTYPESH
#define HDFTYPESH

#include <hdf5.h>
#include <serror.h>

namespace skn
{

template <class T> hid_t hdf_mem_type() {serror("Type not implemented");return hid_t();}
template <> hid_t hdf_mem_type<double>() { return H5T_NATIVE_DOUBLE; }
template <> hid_t hdf_mem_type<float>() { return H5T_NATIVE_FLOAT; }
template <> hid_t hdf_mem_type<int>() { return H5T_NATIVE_INT; }

template <class T> hid_t hdf_file_type() {serror("Type not implemented");return hid_t();}
template <> hid_t hdf_file_type<double>() { return H5T_IEEE_F64LE; }
template <> hid_t hdf_file_type<float>() { return H5T_IEEE_F32LE; }
template <> hid_t hdf_file_type<int>() { return H5T_STD_I32LE; }

hid_t dtype2hdf(const std::string & dtype)
{
	if(dtype == "double") return hdf_file_type<double>();
	if(dtype == "float")  return hdf_file_type<float>();
	if(dtype == "int")    return hdf_file_type<int>();
	serror("Unknown data type %s", dtype.c_str());
	return hdf_file_type<double>();
}

}

#endif
