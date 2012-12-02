/* This file is part of libtrevisan, a modular implementation of
   Trevisan's randomness extraction construction.

   Copyright (C) 2011-2012, Wolfgang Mauerer <wm@linux-kernel.net>

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with libtrevisan. If not, see <http://www.gnu.org/licenses/>. */

#include"timing.h"
#include<iostream>

#define NS_PER_MS (1000*1000)
#define MS_PER_SEC 1000

static double conv_factor;

void init_timekeeping() {
#ifdef __MACH__
	mach_timebase_info_data_t timebase;
	mach_timebase_info(&timebase);
	conv_factor = (double)timebase.numer / (double)timebase.denom;
#endif
}

int timestamp_subtract(meas_t *result, meas_t *t2, meas_t *t1) {
#ifdef __MACH__
	*result = *t2 - *t1;
	return (*result < 0);
#else
	long long diff = (t2->tv_nsec + NSEC_PER_SEC*t2->tv_sec) - 
		(t1->tv_nsec + NSEC_PER_SEC*t1->tv_sec);

	result->tv_sec = diff/NSEC_PER_SEC;
	result->tv_nsec = diff % NSEC_PER_SEC;
	
	if (diff < 0) {
		std::cout << "WARNING: Negative time diff" << std::endl;
		std::cout << "(did you use the high-precision mechanism" << std::endl;
		std::cout << " over extended time periods?)" << std::endl;
	}

	return (diff<0);
#endif
}

void measure(meas_t *timestamp) {
#ifdef __MACH__
	*timestamp = mach_absolute_time();
#else
	clock_gettime(CLOCK_MONOTONIC, timestamp);
#endif
}

long delta_to_ms(meas_t delta) {
#ifdef __MACH__
	double delta_ns = delta*conv_factor;
	return (long)delta_ns/NS_PER_MS;
#else
	return delta.tv_sec*MS_PER_SEC + delta.tv_nsec/NS_PER_MS;
#endif
}

long delta_to_ns(meas_t delta) {
#ifdef __MACH__
	double delta_ns = delta*conv_factor;
	return delta_ns;
#else
	return delta.tv_sec*NSEC_PER_SEC + delta.tv_nsec;
#endif
}

float delta_to_s(meas_t delta) {
#ifdef __MACH__
	float delta_ns = delta*conv_factor;
	return delta_ns/NSEC_PER_SEC;
#else
	return (float)delta.tv_sec + delta.tv_nsec/(float)NSEC_PER_SEC;
#endif
}
