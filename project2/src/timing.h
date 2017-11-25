#ifndef TIMING_H
#define TIMING_H

#ifdef _WIN32
	#include <windows.h>
	#include <winbase.h>
#endif

#if defined(__linux__) || defined(__APPLE__)
	#include <sys/time.h>
#endif

// Timer() function, returns elapsed 'wall-clock' time in seconds
// Time zero is set the first time Timer() is called.
double Timer()
{
#ifdef _WIN32
	static LARGE_INTEGER freq, li0;
    static bool dummy1 = QueryPerformanceFrequency(&freq);
    static bool dummy2 = QueryPerformanceCounter(&li0);
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return static_cast<double>(li.QuadPart-li0.QuadPart)/static_cast<double>(freq.QuadPart);
#endif

#if defined(__linux__) || defined(__APPLE__)
	static timeval tvStart;
    static int dummy = gettimeofday(&tvStart, 0);
	timeval tv;
	gettimeofday(&tv, 0);
	return double(tv.tv_sec-tvStart.tv_sec)+double(tv.tv_usec-tvStart.tv_usec)*0.000001;
#endif
}

#endif
