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

#ifndef DEBUG_H
#define DEBUG_H

// TODO: Stolen from wpa_supplicant
#ifdef __GNUC__
#define PRINT_FORMAT(a,b) __attribute__ ((format (printf, (a), (b))))
#else
#define PRINT_FORMAT(a,b)
#endif

void debug_msg(int level, const char *fmt, ...);

#include "debug_levels.h"

#endif
