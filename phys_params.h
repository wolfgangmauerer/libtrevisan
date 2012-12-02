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

#ifndef PHYS_PARAMS_H
#define PHYS_PARAMS_H

// All physical parameters if the Trevisan process are collected here.
struct phys_params {
	uint64_t n; // Input length (bits)
	uint64_t m; // Output length (bits)
	
	double alpha; // Source entropy ratio
	double eps;   // 1-bit error probability
	
	// Parameters for specific extractors (see the paper for details)
	double lu_nu;
};

#endif
