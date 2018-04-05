/********************************************************************************************
* Hardware-based random number generation function using /dev/urandom
* This file is copied from the random.c file inside SIKE library developed by Microsoft Research
*********************************************************************************************/

#include "random.h"
#include <stdlib.h>
#if defined(__WINDOWS__)
#include <windows.h>
#include <bcrypt.h>
#elif defined(__LINUX__)
#include <unistd.h>
#include <fcntl.h>
static int lock = -1;
#endif


static __inline void delay(unsigned int count)
{
	while (count--) {}
}
#define passed 0 
#define failed 1

int randombytes(unsigned char* random_array, unsigned long long nbytes)
{ // Generation of "nbytes" of random values

#if defined(__WINDOWS__)   
	if (!BCRYPT_SUCCESS(BCryptGenRandom(NULL, random_array, (unsigned long)nbytes, BCRYPT_USE_SYSTEM_PREFERRED_RNG))) {
		return failed;
	}

#elif defined(__LINUX__)
	int r, n = (int)nbytes, count = 0;

	if (lock == -1) {
		do {
			lock = open("/dev/urandom", O_RDONLY);
			if (lock == -1) {
				delay(0xFFFFF);
			}
		} while (lock == -1);
	}

	while (n > 0) {
		do {
			r = read(lock, random_array + count, n);
			if (r == -1) {
				delay(0xFFFF);
			}
		} while (r == -1);
		count += r;
		n -= r;
	}
#endif

	return passed;
}