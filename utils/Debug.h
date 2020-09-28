#ifndef MV_Debug_H_

#ifdef NDEBUG
#define MV_ASSERT(expression)
#define MV_LOGD(x)
#else
#include <iostream>
#include <cassert>

#define MV_ASSERT(expression) assert(expression);
#define MV_CLOGD(x) std::clog << "Debug: " << (x) << std::endl;
#define MV_LOGD(...) printf(__VA_ARGS__); printf("\n");	
#endif

#endif // MV_Debug_H_