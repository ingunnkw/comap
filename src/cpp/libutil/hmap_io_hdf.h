#ifndef HMAPIOHDFH
#define HMAPIOHDFH

/* This header provides read_hdf_map and write_hdf_map */

#include <hdf_types.h>
#include <hmap.h>

namespace skn {

template <class T>
void write_scalar(hid_t f, const std::string & name, const T & data)
{
	hid_t space = H5Screate_simple(0, NULL, NULL);
	hid_t set   = H5Dcreate(f, name.c_str(), hdf_file_type<T>(), space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(set, hdf_mem_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
}
template <class T>
void read_scalar(hid_t f, const std::string & name, T & data)
{
	hid_t set = H5Dopen(f, name.c_str(), H5P_DEFAULT);
	hid_t space = H5Dget_space(set);
	if(H5Sget_simple_extent_ndims(space) != 0) serror("Set is not a scalar!");
	H5Dread(set, hdf_mem_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
}

template <typename T>
HMap<T> read_hmap_hdf(const std::string & fname)
{
	hid_t file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	// Get general properties
	int nside, ordering, fullsky;
	read_scalar(file, "nside", nside);
	read_scalar(file, "ordering", ordering);
	read_scalar(file, "fullsky", fullsky);

	// Get the extents
	hid_t mapset = H5Dopen(file, "maps", H5P_DEFAULT);
	hid_t mapspace = H5Dget_space(mapset);
	Vector<int,3> ext;
	Vector<hsize_t,3> hext;
	if(H5Sget_simple_extent_ndims(mapspace) != 3)
		serror("'maps' array does not have 3 dimensions!");
	H5Sget_simple_extent_dims(mapspace, &hext[0], NULL);
	for(int i = 0; i < 3; i++) ext[i] = hext[i];

	// Get the maps
	Grid<T,3> data(ext);
	H5Dread(mapset, hdf_mem_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

	// And the pixels
	if(!fullsky)
	{
		hid_t pixset = H5Dopen(file, "pixels", H5P_DEFAULT);
		hid_t pixspace = H5Dget_space(pixset);
		hsize_t len;
		if(H5Sget_simple_extent_ndims(pixspace) != 1)
			serror("'pixels' array does not have 1 dimension!");
		H5Sget_simple_extent_dims(pixspace, &len, NULL);
		if(len != ext(2))
			serror("Inconsistent npix between pixels an maps arrays! (%d vs %d)",
				len, ext(2));
		Grid<int,1> pixels(len);
		H5Dread(pixset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pixels[0]);

		return HMap<T>(nside, ordering, pixels, data);
	}
	else return HMap<T>(nside, ordering, data);
}

template <typename T>
void write_hmap_hdf(const std::string & fname, const HMap<T> & res, hid_t pred)
{
	hid_t file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	write_scalar(file, "nside", res.get_nside());
	write_scalar(file, "ordering", res.get_ordering());
	write_scalar(file, "fullsky", (int)res.is_fullsky());

	hsize_t ext[3]; for(int i = 0; i < 3; i++) ext[i] = res.size(i);
	hid_t mapspace = H5Screate_simple(3, ext, NULL);
	hid_t mapset   = H5Dcreate(file, "maps", pred, mapspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(mapset, hdf_mem_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &res[0]);

	if(!res.is_fullsky())
	{
		hsize_t len = res.size(2);
		hid_t pixspace = H5Screate_simple(1, &len, NULL);
		hid_t pixset   = H5Dcreate(file, "pixels", H5T_STD_I32LE, pixspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(pixset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &res.get_pixels()[0]);
	}
}

template <typename T>
void write_hmap_hdf(const std::string & fname, const HMap<T> & m, const char * dtype = 0)
{
	if(!dtype) write_hmap_hdf(fname, m, hdf_file_type<T>());
	else write_hmap_hdf(fname, m, dtype2hdf(dtype));
}

}

#endif
