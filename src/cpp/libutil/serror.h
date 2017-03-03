#ifndef MYERRORFUN
#define MYERRORFUN

#include <cstdarg>
#include <string>
namespace skn
{
	struct Error { Error(const std::string & m):msg(m) {} std::string msg; };
	void serror(const char * fmt, ...);
}
#endif
