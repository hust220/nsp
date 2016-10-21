#pragma once

#define FREE_ARR(a) do { \
	if (a != NULL) delete [] a; \
	a = NULL; \
} while (0)

#define FREE_OBJ(a) do { \
	if (a != NULL) delete a; \
	a = NULL; \
} while (0)

namespace jian {
	using uint = unsigned int;
}
