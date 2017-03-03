#include <serror.h>
#include <cstdio>
#include <cstdarg>
using namespace std;
namespace skn {
	void serror(const char * fmt, ...)
	{
		using namespace std;
		char buffer[0x1000];
		va_list ap;
		va_start(ap, fmt);
		vsprintf(buffer, fmt, ap);
		va_end(ap);
		throw Error(buffer);
	}
}
