#include <stimer.h>
#include <sys/time.h>

namespace skn {
	double wall_time() { timeval tv; gettimeofday(&tv,0); return tv.tv_sec + 1e-6*tv.tv_usec; }
}
