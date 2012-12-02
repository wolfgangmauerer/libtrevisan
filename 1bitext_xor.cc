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

// A 1-bit extractor based on XOR
// Algorithmic description and parameter calculations by Ch. Portmann

#include<iostream>
#include<fstream>
#include<cstddef>
#include<stdlib.h>
#include<RInside.h>
#include "timing.h"
#include "utils.hpp"
#include "1bitext_xor.h"

using namespace std;

vertex_t bitext_xor::num_random_bits() {
	if (l ==  std::numeric_limits<vertex_t>::max()) {
		cerr << "(XOR extractor) Internal error: num_random_bits() called for "
		     << "without valid value for parameter l" << endl;
		exit(-1);
	}
	// The number of random bits per run is l*log(n) (see the comment in
	// extract why we compute with nb(n-1))
	return l*numbits<vertex_t>(pp.n-1);
}

void bitext_xor::compute_l() {
	SEXP ans;
	stringstream call;

	call << "do.opt.xor(" << pp.alpha << ", " << mu << ", " << r << ", "
	     << pp.n << ", " << pp.eps << ")";

	r_interp->parse_eval(call.str(), ans);
	l = Rcpp::as<vertex_t>(ans);
}

uint64_t bitext_xor::compute_k() {
	SEXP ans;
	stringstream call;
	call << "gamma.func(" << pp.alpha << ", " << mu << ", " << r << ", "
	     << pp.n << ", " << pp.eps << ")" << endl;

	r_interp->parse_eval(call.str(), ans);
	long double gamma = Rcpp::as<long double>(ans);

	long double eps_prime = sqrt(pp.eps)/(1+sqrt(2));
	vertex_t k = ceil(gamma*pp.n + r*pp.m + 6*log2((1+sqrt(2))/eps_prime) +
			  log2(4.0/3.0));

	return(k);
}

bool bitext_xor::extract(void *initial_rand) {
	unsigned short bits_per_val = numbits<uint64_t>(pp.n-1);
	unsigned short bits_per_datum = BITS_PER_TYPE(datum_t);
	datum_t idx;
	bool res = 0;
	local_bf.set_raw_data(initial_rand, num_random_bits());

#ifdef EXPENSIVE_SANITY_CHECKS
	if (l ==  std::numeric_limits<vertex_t>::max()) {
		cerr << "(XOR extractor) Internal error: extract() called for "
		     << "without valid value for parameter l" << endl;
		exit(-1);
	}
#endif

	// NOTE: Since the number of bits representable in a machine with
	// word length w is 8*2^w= 2^(w+3), we can simplify the bit selection
	// if sizeof(elem_t) \geq w+3: It is guaranteed that all bits fit
	// into one instance of the type in principle.

	// Technically, this amounts to the fact that elem_t must provide
	// enough bits to hold a complete index, that is, n bits.
	if (bits_per_datum < bits_per_val) {
		fprintf(stderr, "Internal error: Assumption "
			"bits_per_datum < bits_per_val failed!\n");
		exit(-1);
	}

	// Iterate along the given initial randomness, take nb(n-1)=bits_per_val bits in
	// each step, and select the bit in n with this index. Compute the
	// cumulative XOR value over all bits.
	for (uint64_t i = 0; i < l; i++) {
		local_bf.get_bit_range(i*bits_per_val, (i+1)*bits_per_val-1, &idx);

		// Take the bit indexed by idx of the global randomness
		// and xor it to the result
		res ^= b.get_bit(idx);
	}

	return res;
}
