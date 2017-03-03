#ifndef HPIXDATA
#define HPIXDATA

#include <string>
#include <vector>
#include <arr.h>
#include <stdio.h>

namespace skn {

std::string trypath(const std::vector<std::string> & paths, const std::string & filename);
std::vector<std::string> hpix_paths();
arr<double> get_hpix_array(int nside, char * fmt, double add = 0);
const arr<double> & get_weights(int nside, int lmax = 0);
const arr<double> & get_window(int nside);

}

#endif
