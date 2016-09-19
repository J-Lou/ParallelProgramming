///////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2005, 2006, 2007, 2008, 2009
// University of Rochester
// Department of Computer Science
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of the University of Rochester nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef __HRTIME_H__
#define __HRTIME_H__

#if defined(__i386__) && !defined(__APPLE__) && !defined(_MSC_VER)
// gethrtime implementation by Kai Shen for x86 Linux

#ifdef __linux__
#include <stdio.h>
#include <string.h>
#include <assert.h>

// get the number of CPU cycles per microsecond from Linux /proc filesystem
// return < 0 on error
inline double getMHZ_x86(void)
{
    double mhz = -1;
    char line[1024], *s, search_str[] = "cpu MHz";
    FILE* fp;

    // open proc/cpuinfo
    if ((fp = fopen("/proc/cpuinfo", "r")) == NULL)
        return -1;

    // ignore all lines until we reach MHz information
    while (fgets(line, 1024, fp) != NULL) {
        if (strstr(line, search_str) != NULL) {
            // ignore all characters in line up to :
            for (s = line; *s && (*s != ':'); ++s);
            // get MHz number
            if (*s && (sscanf(s+1, "%lf", &mhz) == 1))
                break;
        }
    }

    if (fp != NULL)
        fclose(fp);
    return mhz;
}
#endif // __linux__

#if defined(__FreeBSD__) || defined(__OpenBSD__)
#include <sys/param.h>
#include <sys/sysctl.h>

inline double getMHZ_x86()
{
   int mhz = -1;
   size_t mhz_sz = sizeof(mhz);
#if defined(__FreeBSD__)
   int cmd[2];
   size_t two = 2;
   sysctlnametomib("hw.clockrate", cmd, &two);
#elif defined(__OpenBSD__)
   int cmd[] = {CTL_HW, HW_CPUSPEED};
#endif
   sysctl(cmd, sizeof(cmd)/sizeof(int), &mhz, &mhz_sz, NULL, 0);
   return (double)mhz;
}
#endif // FreeBSD || OpenBSD

// get the number of CPU cycles since startup using rdtsc instruction
inline unsigned long long gethrcycle_x86()
{
    unsigned int tmp[2];

    asm ("rdtsc"
         : "=a" (tmp[1]), "=d" (tmp[0])
         : "c" (0x10) );
    return (((unsigned long long)tmp[0] << 32 | tmp[1]));
}

// get the elapsed time (in nanoseconds) since startup
inline unsigned long long getElapsedTime()
{
    static double CPU_MHZ = 0;
    if (CPU_MHZ == 0)
        CPU_MHZ = getMHZ_x86();
    return (unsigned long long)(gethrcycle_x86() * 1000 / CPU_MHZ);
}
#endif // i386

#if defined(__linux__) && defined(__ia64__)

// get the number of CPU cycles per microsecond from Linux /proc filesystem
// return < 0 on error
inline double getMHZ_ia64(void)
{
    double mhz = -1;
    char line[1024], *s, search_str[] = "itc MHz";
    FILE* fp;

    // open proc/cpuinfo
    if ((fp = fopen("/proc/cpuinfo", "r")) == NULL)
        return -1;

    // ignore all lines until we reach MHz information
    while (fgets(line, 1024, fp) != NULL) {
        if (strstr(line, search_str) != NULL) {
            // ignore all characters in line up to :
            for (s = line; *s && (*s != ':'); ++s);
            // get MHz number
            if (*s && (sscanf(s+1, "%lf", &mhz) == 1))
                break;
        }
    }

    if (fp != NULL)
        fclose(fp);
    return mhz;
}

// get cpu cycles by reading the itc
inline unsigned long gethrcycle_ia64()
{
    unsigned long tmp;
    asm volatile("mov %0=ar.itc;;":"=r"(tmp));
    return tmp;
}

// get the elapsed time (in nanoseconds) since startup
//
// We're returning a ull so that we don't have to modify the call sites, but
// really on ia64 a ul would be sufficient.
inline unsigned long long getElapsedTime()
{
    static double CPU_MHZ = 0;
    if (CPU_MHZ == 0)
        CPU_MHZ = getMHZ_ia64();
    return (unsigned long long)(gethrcycle_ia64() * 1000 / CPU_MHZ);
}
#endif

#ifdef __sparc__

// use gethrtime() on SPARC/SOLARIS
inline unsigned long long getElapsedTime()
{
    return gethrtime();
}

#endif // SPARC

#ifdef __APPLE__
// Darwin timer

// This code is based on code available through the Apple developer's
// connection website at http://developer.apple.com/qa/qa2004/qa1398.html

#include <mach/mach_time.h>

inline unsigned long long getElapsedTime()
{
    static mach_timebase_info_data_t sTimebaseInfo;

    // If this is the first time we've run, get the timebase.
    // We can use denom == 0 to indicate that sTimebaseInfo is
    // uninitialised because it makes no sense to have a zero
    // denominator as a fraction.

    if ( sTimebaseInfo.denom == 0 )
        (void)mach_timebase_info(&sTimebaseInfo);

    return mach_absolute_time() * sTimebaseInfo.numer / sTimebaseInfo.denom;
}

#endif

#ifdef _MSC_VER
#include <windows.h>
#pragma comment(lib, "winmm.lib")
inline unsigned long long getElapsedTime()
{
    unsigned long long ret;
    QueryPerformanceCounter((LARGE_INTEGER*)&ret);
    return ret;
    // return timeGetTime() * 1000;
}

inline int gettimeofday(struct timeval* tp, void* tzp)
{
    DWORD t;
    t = timeGetTime();
    tp->tv_sec = t / 1000;
    tp->tv_usec = t % 1000;
    // 0 indicates success
    return 0;
}
#endif

#ifdef _AIX
// On AIX use read_real_time()
inline unsigned long long getElapsedTime()
{

  timebasestruct_t r;

  // First we read the real time
  read_real_time(&r, sizeof(r));

  // The following line converts to POWER format if running on a
  // PPC601.  The manpage recommends using it at any rate, probably
  // since newer hardware could conceivably change the format.
  time_base_to_time(&r, sizeof(r));

  // Now, r.tb_high is seconds, and r.tb_low is nanoseconds
  unsigned long long sec = r.tb_high, usec = r.tb_low;

  return sec * 1000000000ULL + usec;

}
#endif

// -- sleep_ms ---------------------------------------------------------- //
//  Sleep is different on different OSes.
#ifdef _MSC_VER
inline void sleep_ms(unsigned int ms) { Sleep(ms); }
#else
#include <unistd.h>
inline void sleep_ms(unsigned int ms) { usleep(ms*1000); }
#endif

// -- yield ---------------------------------------------------------- //
//  Yield is different on different OSes.
#if defined(__sparc__)
inline void yield_cpu() { yield(); }
#elif defined(__APPLE__)
inline void yield_cpu() { sched_yield(); }
#elif defined(_MSC_VER)
#include <Winbase.h>
inline void yield_cpu() { SwitchToThread(); }
#elif defined(_AIX)
#include <sched.h>
inline void yield_cpu() { sched_yield(); }
#else
#include <pthread.h>
inline void yield_cpu() { pthread_yield(); }
#endif

#endif // __HRTIME_H__
