/***************************************************************
   File:     papi_timer.c

   Function: get_cur_time()

   Return:   the time since some arbitrary starting point in seconds
             as a double value.
             PAPI timers use the most accurate timers available 
             on the platform in use.

   Note:     must be used on the platform where PAPI is installed 
             (e.g. torc1), and remember to call PAPI_library_init() 
             first.

****************************************************************/

#include <papi.h>

double get_cur_time() {
  return PAPI_get_real_usec() / 1000000.0;
} 


