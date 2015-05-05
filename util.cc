// Author: Mingcheng Chen (linyufly@gmail.com)

#include "util.h"

#include <cstdio>
#include <cstdarg>
#include <cstdlib>

#ifdef IN_VC

#include <time.h>
#include <windows.h>

#else

#include <sys/time.h>

#endif

namespace {

void report_va_error(const char *format, va_list args) {
  vfprintf(stderr, format, args);

  exit(EXIT_FAILURE);
}

#ifdef IN_VC

// Copy from http://social.msdn.microsoft.com/Forums/vstudio/en-US/
// 430449b3-f6dd-4e18-84de-eebd26a8d668/gettimeofday?forum=vcgeneral
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag;
 
  if (NULL != tv)
  {
    GetSystemTimeAsFileTime(&ft);
 
    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;
 
    /*converting file time to unix epoch*/
    tmpres -= DELTA_EPOCH_IN_MICROSECS; 
    tmpres /= 10;  /*convert into microseconds*/
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }
 /*
  if (NULL != tz)
  {
    if (!tzflag)
    {
      _tzset();
      tzflag++;
    }
    tz->tz_minuteswest = _timezone / 60;
    tz->tz_dsttime = _daylight;
  }
 */
  return 0;
}

#endif

}

void report_error(const char *format, ...) {
  va_list args;
  va_start(args, format);
  report_va_error(format, args);
  va_end(args);
}

double get_time_in_seconds() {
  timeval currTime;
  gettimeofday(&currTime, 0);
  return currTime.tv_sec + currTime.tv_usec * 1e-6;
}

