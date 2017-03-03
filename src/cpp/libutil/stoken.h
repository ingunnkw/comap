#ifndef STOKKYO
#define STOKKYO

#include <string>
#include <vector>

namespace skn
{

std::vector<std::string> tokenize(const std::string & line, const std::string & sep);
std::vector<std::string> tokenize_full(const std::string & line, const std::string & sep);

}

#endif
