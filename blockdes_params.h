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

#ifndef BLOCKDES_PARAMS_H
#define BLOCKDES_PARAMS_H

#include <RInside.h>
#include <vector>
#include <cstdint>
#include "R_interp.h"

typedef std::vector<std::vector<uint64_t> > blockdes_properties;

// Interface to R for computing block design parameters
class blockdes_params {
public:
	blockdes_params(class R_interp *_r_interp) : r_interp(_r_interp) {
		init_functions();		
	};

	Rcpp::IntegerMatrix compute_blocks(double r, uint64_t m, uint64_t t);

private:
	R_interp *r_interp;
	
	void init_functions();
};

#endif
