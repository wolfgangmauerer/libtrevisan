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

#ifndef TIMING_H
#define TIMING_H

#include <time.h>
#include <string>

#ifdef __MACH__
#include <mach/mach_time.h>
#include <tr1/cstdint>
typedef uint64_t meas_t;
#else
typedef struct timespec meas_t;
#define NSEC_PER_SEC    1000000000L
#endif

void init_timekeeping();
int timestamp_subtract(meas_t *result, meas_t *x, meas_t *y);
void measure(meas_t *timestamp);
long delta_to_ms(meas_t delta);
long delta_to_ns(meas_t delta);
float delta_to_s(meas_t delta);

#endif
