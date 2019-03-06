#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

/* This file contatins the following timing functions:

    rdtscp(i)
    : Read Time-Stap Counter and Processor ID

    timer_start(timer_id)
    timer_end(timer_id, sec, nsec)
    : Uses the RTC of the system to measure time

    cpseed(polls)
    : Uses rdtscp and timer_start and timer_end to caclulate the
    : speed of the system's cpu in hertz.

*/

/* In Fortran, use the following as an interface for rdtscp:
    use iso_c_binding, only : c_long
 
    interface
        subroutine rdtscp(i) bind(C)
           use iso_c_binding, only : c_long
           integer (c_long), intent(out) :: i
        end subroutine rdtscp
    end interface
 
    integer (c_long) :: tsc_start, tsc_end 
*/
    
/* In Fortran, use the following as an interfaces for rtc timer:
    use iso_c_binding, only : c_int
 
    interface
        subroutine timer_start(timer_id) bind(C)
           use iso_c_binding, only : c_int
           integer (c_int), intent(in), value :: timer_id
        end subroutine timer_start

        subroutine timer_stop(timer_id, sec, nsec) bind(C)
           use iso_c_binding, only : c_int
           integer (c_int), intent(in), value :: timer_id
           integer (c_int), intent(out) :: sec, nsec
        end subroutine timer_stop
    end interface
 
    integer (c_int) :: timer_id, sec, nsec 
*/

/* In Fortran, use the following as an interface for cspeed:

   use iso_c_binding, only : c_float, c_int
    
    interface
        real(c_float) function cspeed(polls) BIND(C)
            use iso_c_binding, only : c_float, c_int
            integer(c_int) :: polls
        end function cspeed
    end interface

    integer (c_int) :: polls
    real (c_float) :: clockSpeed
*/

void rdtscp( long *i )
{
	unsigned rax, rdx;
	asm volatile ("RDTSCP\n\t"
			"mov %%edx, %0\n\t"
			"mov %%eax, %1\n\t": "=r" (rdx), "=r" (rax));
	*i = ((unsigned long)rdx << 32) + rax;
}

#define MAX_TIMERS 10

#ifdef GETTIMEOFDAY
#include <sys/time.h>
#endif

#ifdef __MACH__
#include <mach/mach.h>
#include <mach/mach_time.h>
#include <unistd.h>
#endif

#ifdef __linux__
#include <time.h>
#endif

#ifdef GETTIMEOFDAY
struct timeval start_time[MAX_TIMERS];
struct timeval end_time[MAX_TIMERS];
#endif

#ifdef __MACH__
uint64_t start_time[MAX_TIMERS];
uint64_t end_time[MAX_TIMERS];
#endif

#ifdef AIX
timebasestruct_t start_time[MAX_TIMERS];
timebasestruct_t end_time[MAX_TIMERS];
#endif

#ifdef __linux__
struct timespec start_time[MAX_TIMERS];
struct timespec end_time[MAX_TIMERS];
#endif

void timer_start(int n)
{
#ifdef GETTIMEOFDAY
   gettimeofday(&start_time[n], NULL);
#endif

#ifdef __MACH__
   start_time[n] = mach_absolute_time();
#endif

#ifdef AIX
   read_real_time(&start_time[n], TIMEBASE_SZ);
#endif

#ifdef __linux__
   clock_gettime(CLOCK_MONOTONIC_RAW, &start_time[n]);
#endif
}

void timer_stop(int n, int *secs, int *n_secs)
{
#ifdef GETTIMEOFDAY
   gettimeofday(&end_time[n], NULL);
  
   *secs   = (int)(end_time[n].tv_sec - start_time[n].tv_sec);
   *n_secs = (int)(end_time[n].tv_usec - start_time[n].tv_usec) * 1000;

   if (*n_secs < 0)  {
      *secs   -= 1;
      *n_secs += 1000000000;
   }
#endif

#ifdef __MACH__
   uint64_t elapsed, elapsedNano;
   static mach_timebase_info_data_t sTimebaseInfo;

   end_time[n] = mach_absolute_time();

   elapsed = end_time[n] - start_time[n];

    if ( sTimebaseInfo.denom == 0 ) {
        (void) mach_timebase_info(&sTimebaseInfo);
    }

    // Do the maths. We hope that the multiplication doesn't 
    // overflow; the price you pay for working in fixed point.

    elapsedNano = elapsed * sTimebaseInfo.numer / sTimebaseInfo.denom;


   *secs   = (int)(elapsedNano / 1000000000);
   *n_secs = (int)(elapsedNano % 1000000000);
#endif

#ifdef AIX
   read_real_time(&end_time[n], TIMEBASE_SZ);
   time_base_to_time(&start_time[n], TIMEBASE_SZ);
   time_base_to_time(&end_time[n], TIMEBASE_SZ);

   *secs = end_time[n].tb_high - start_time[n].tb_high;
   *n_secs = end_time[n].tb_low - start_time[n].tb_low;

   if (*n_secs < 0)  {
      *secs   -= 1;
      *n_secs += 1000000000;
   }
#endif

#ifdef __linux__
   clock_gettime(CLOCK_MONOTONIC_RAW, &end_time[n]);

   *secs = (int)(end_time[n].tv_sec - start_time[n].tv_sec);
   *n_secs = (int)(end_time[n].tv_nsec - start_time[n].tv_nsec);

   if (*n_secs < 0)  {
      *secs   -= 1;
      *n_secs += 1000000000;
   }
#endif
}

/* cspeed - Determine the machines CPU Speed
 * Calculate the CPU speed of the machine by dividing the number of
 * clock cycles between a time frame
 *
 * This calculation is only realiable on intel X86 CPU's that contain the
 * flag: constant_tsc. This flag can be checked by viewing the contents of
 * `/proc/cpuinfo`.
 *
 */
float cspeed(int *polls){

    long clockSpeed = 0;
    long tsc_start = 0;
    long tsc_end = 0;
    long tsc_elapsed = 0;
    int timer_id = 0;
    int sec = 0;
    int nsec = 0;
    int i;

    for (i = 0; i < *polls; i++){
        rdtscp(&tsc_start);
        timer_start(0);
        sleep(3);
        timer_stop(timer_id, &sec, &nsec);
        rdtscp(&tsc_end);
        tsc_elapsed = tsc_end - tsc_start;
        clockSpeed += ( tsc_elapsed / ( sec + nsec*10E-9 ));
    }

    clockSpeed = clockSpeed / ((float) *polls);

    return roundf(clockSpeed * 10E-9) / 10E-9;
    //return clockSpeed;
}

