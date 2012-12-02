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

#include <RInside.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "blockdes_params.h"

using namespace std;

void blockdes_params::init_functions() {
	string code =
#include "generated/bd_r_embedd.inc"
	r_interp->eval_code(code);
}

Rcpp::IntegerMatrix blockdes_params::compute_blocks(double r, uint64_t m,
						    uint64_t t) {
	SEXP ans;
	stringstream call;
	call << "block.param.as.matrix(" << r << ", " << m << ", " << t << ")";

	r_interp->parse_eval(call.str(), ans); 
	Rcpp::IntegerMatrix M(ans);

	return (M);
}
