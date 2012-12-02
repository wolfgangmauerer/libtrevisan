#include <varargs.h>
#include "debug.h"

extern short debug_level;

void debug_msg(int level, const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	if (level >= debug_level) {
		vprintf(fmt, ap);
		printf("\n");
	}
	va_end(ap);
}
