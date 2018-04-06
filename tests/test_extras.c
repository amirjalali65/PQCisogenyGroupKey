/********************************************************************************************
* Supersingular Isogeny Group Key Agreement Library
*
* Abstract: utility functions for testing and benchmarking
*********************************************************************************************/

#include "test_extras.h"
#if (OS_TARGET == OS_LINUX) && (TARGET == TARGET_ARM || TARGET == TARGET_ARM64)
    #include <time.h>
#endif
#include <stdlib.h>


int64_t cpucycles(void)
{  //Access system counter for benchmarking
    unsigned int hi, lo;

    asm volatile ("rdtsc\n\t" : "=a" (lo), "=d"(hi));
    return ((int64_t)lo) | (((int64_t)hi) << 32);
#if (OS_TARGET == OS_LINUX) && (TARGET == TARGET_ARM || TARGET == TARGET_ARM64)
    struct timespec time;

    clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec*1e9 + time.tv_nsec);
#else
    return 0;            
#endif
}
