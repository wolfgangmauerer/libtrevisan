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

// A one-bit extractor based on g-regular graphs as described by Lu
#include<iostream>
#include<fstream>
#include<cstddef>
#include<cstdlib>
#include<gmp.h>
#include "timing.h"
#include "utils.hpp"
#include "1bitext_expander.h"

#ifndef NO_DEBUG
  #include "debug.h"
extern int debug_level;
#else
#include "debug_levels.h"
int debug_level = 0;
void debug_msg(int level, const char *fmt, ...) { }
#endif

using namespace std;

void bitext_expander::infer_params() {
	SEXP ans;
	stringstream call;

	call << "do.compute.lu(" << nu << ", " << pp.m << ", " << pp.eps << ", "
	     << lambda0 << ")";
	r_interp->parse_eval(call.str(), ans);
	Rcpp::DataFrame res(ans);
	w = Rcpp::as<long double>(res(0));
	c = Rcpp::as<vertex_t>(res(1));
	l = Rcpp::as<vertex_t>(res(2));
}

uint64_t bitext_expander::compute_k() {
	return(h(nu)*pp.n + r*pp.m + 6.0*log((2.0+sqrt(2.0))/pp.eps) - 2.0);
}

void bitext_expander::set_input_data(void *global_rand, struct phys_params &pp) {
	bitext::set_input_data(global_rand, pp);
	nu = pp.lu_nu;
	infer_params();

	// Ensure that n is of the form b^2 for some b \in
	// \mathbb{N}^{+}, b even (b is sqrt_n in the following)
	mpz_t n_gmp; 
	mpz_init_set_ui(n_gmp, pp.n);

	mpz_t sqrt_n_gmp;
	mpz_init(sqrt_n_gmp);

	if (mpz_root(sqrt_n_gmp, n_gmp, 2) == 0) {
		cerr << "(Lu bit extractor) Internal error: n != b^2" << endl;
		cerr << "(n=" << pp.n << ")" << endl;
		exit(-1);
	}

	if (mpz_even_p(sqrt_n_gmp) == 0) {
		// NOTE: We could also extend the scheme to odd
		// bit numbers, but that would unecessarily
		// complicate things in the graph calculations
		// without benefit.
		cerr << "(Lu bit extractor) Error: Only even b (for n=b^2)is supported!"
		     << endl;
		exit(-1);
	}

	if (mpz_fits_ulong_p (sqrt_n_gmp) == 0) {
		cerr << "(Lu bit extractor) Internal error: m does not fit into vertex_t!"
		     << endl;
		exit(-1);
	}

	sqrt_n = mpz_get_ui(sqrt_n_gmp);

	// m needs half the amount of bits of n because n = m^2
	// => m = sqrt(n) => log(m) = log(n^{1/2}) = 1/2*log(n)
	sqrt_n_bits = numbits<vertex_t>(sqrt_n-1); // sqrt_n values fit into [0,sqrt_n-1]
	n_bits = 2*sqrt_n_bits;
	zb_zb_mask = (static_cast<vertex_t>(1) << sqrt_n_bits) - 1;

	if (debug_level >= INFO) {
		cerr << "(Lu extractor) bits(n): " << n_bits << ", bits(sqrt(n)): "
		     << sqrt_n_bits << endl;
	}

	// Due diligence
	mpz_clear(n_gmp);
	mpz_clear(sqrt_n_gmp);

	// Compute offsets for the different portions of the initial randomness
	// (the actual pointers can differ, but the offsets are invariant)
	additional_randomness = 0;

	// Offsets are computed in terms of edge_datum_t, and if the amount
	// of bits required for one component is not evenly divisible
	// by sizeof(edge_datum_t), we round up by one instance of edge_datum_t --
	// thus the addition by one.
	if (n_bits % BITS_PER_TYPE(edge_datum_t) != 0) {
		additional_randomness =
			n_bits % BITS_PER_TYPE(edge_datum_t);
		walk_bits_offset = n_bits/BITS_PER_TYPE(edge_datum_t) + 1;
	} else {
		walk_bits_offset = n_bits/BITS_PER_TYPE(edge_datum_t);
	}

	if ((c*(l-1)*bits_per_edge) % BITS_PER_TYPE(edge_datum_t) != 0) {
		additional_randomness +=
			(c*(l-1)*bits_per_edge) % BITS_PER_TYPE(edge_datum_t);
		select_bits_offset = (c*(l-1)*bits_per_edge)/BITS_PER_TYPE(edge_datum_t) + 1;
	} else {
		select_bits_offset = (c*(l-1)*bits_per_edge)/BITS_PER_TYPE(edge_datum_t);
	}
}

inline vertex_t_s bitext_expander::do_mod(vertex_t_s a, vertex_t_s mod) {
	vertex_t_s res = a % mod;

	if (res < 0)
		res += mod;
	
	return res;
}

size_t bitext_expander::multiple_of(size_t num, size_t mult) {
	// Find the smallest n such that n*mult >= num,
	// and return n*mult
	size_t n = num / mult;
	
	if (num % mult != 0)
		return (n+1)*mult;
	
	return n*mult;
}

////////////////////////////////////////////////////////////////////
// NOTE: This function depends on the choice of g-regular graph used
// as basis for the one-bit extractor. Here, we implement the rules
// for a Gabber-Galil expander.
// The number of nodes n=b^2, that is, we operate on 
// \mathbbm{Z}_{b} \times \mathbbm{Z}_{b}.
vertex_t bitext_expander::compute_next_vertex(vertex_t curr_vertex,
					      edge_t next_edge) {
	vertex_t_s x, y;
	// Compute the upper and lower half of Z_b\otimes Z_b
	x = curr_vertex & zb_zb_mask;
	y = (curr_vertex & (zb_zb_mask << sqrt_n_bits)) >> sqrt_n_bits;

	// We do not need to use gmp for the modulo arithmetic here -- 
	// m uses only (at most) half the bits of vertex_t, so we
	// can compute using elementary signed modulo arithmetic
	// TODO: Is this really true for -2INT_MAX % INT_MAX?
	// And likewise for -(INT_MAX+1) % INT_MAX.
	// TODO: This needs to be tested

//	cout << "    Current vertex: " << curr_vertex << ", next edge: "
//		<< next_edge << ", ";
	switch(next_edge) {
	case 0:
		x = do_mod(x + 2*y, sqrt_n);
		break;
	case 1:
		x = do_mod(x - 2*y, sqrt_n);
		break;
	case 2:
		x = do_mod(x + (2*y+1), sqrt_n);
		break;
	case 3:
		x = do_mod(x - (2*y+1), sqrt_n);
		break;
	case 4:
		y = do_mod(y + 2*x, sqrt_n);
		break;
	case 5:
		y = do_mod(y - 2*x, sqrt_n);
		break;
	case 6:
		y = do_mod(y + (2*x+1), sqrt_n);
		break;
	case 7:
		y = do_mod(y - (2*x+1), sqrt_n);
		break;
	default:
		cerr <<  "(Lu bit extractor) Internal error: Edge value=" << next_edge
		     << ", maximal value is 7" << endl;
		exit(-1);
	}

	// Compose the result by concatenating <x,y>
	vertex_t res = x;
	res |=  (y << sqrt_n_bits);
//	cout << "new vertex: " << res << endl;

	return res;
}

// TODO: Check that none of the fixed-precision quantities overflow
// TODO: Include numerical assertions for the invariants

// Determine the required number of random bits for a specific
// parameter set. The required amount may be slightly larger than
// the minimum because of alignment constraints -- this simplifies
// the implementation, but does not cost any significant amount of
// randomness.
vertex_t bitext_expander::num_random_bits() {
	// n_bits aligned by bits(sizeof(edge_datum_t))
	// c*(l-1)*bits_per_edge aligned by bits(sizeof(edge_datum_t))
	// zeta aligned by byte
	vertex_t count;
	
	count = multiple_of(numbits<vertex_t>(pp.n),
			    BITS_PER_TYPE(edge_datum_t));
	count += multiple_of(c*(l-1)*bits_per_edge,
			     BITS_PER_TYPE(edge_datum_t));
	count += multiple_of(l, BITS_PER_BYTE);
	
	if (debug_level >= INFO)
		cout << "Required bits for n=" << pp.n << ", c="
		     << c << ", l=" << l << ": " << count << endl;

	return count;
}

bool bitext_expander::extract(void *initial_rand) {
	unsigned short res = 0; // Output parity: Even -> 0, Odd -> 1
	edge_datum_t *walk_bits;
	void *select_bits;

	// The initial seed is divided into three portions:
	//   - numbits(n) bits to select the initial node (n is the number of vertices)
	//   - c*(l-1)*bits_per_edge bits to perform the random walk (the first of every
	//     c steps is included in the result). Stored in walk_bits.
	//   - l bits to decide if the contribution of the i^{th} node is
	//     included in the result or not. Stored in select_bits
	walk_bits = (edge_datum_t*)initial_rand + walk_bits_offset;
	select_bits = (edge_datum_t*)walk_bits + select_bits_offset;

	if (debug_level >= EXCESSIVE_INFO) {
		cerr << "initial_rand: " << initial_rand << ", walk_bits: " << walk_bits
		     << ", select_bits: " << select_bits << endl;
	}

	if (additional_randomness) {
		if (debug_level >= EXCESSIVE_INFO) {
			cerr << "Information: Requesting " << additional_randomness
			     << " bits more randomness than strictly required because of "
			     << "alignment constraints." << endl;
		}
	}

	if (debug_level >= EXCESSIVE_INFO)
		cout << "Importing " 
		     << multiple_of(l, BITS_PER_BYTE)/BITS_PER_BYTE
		     << " bytes from select_bits" << endl;

	bitfield<uint64_t, uint64_t> select_bits_bf;
	select_bits_bf.set_raw_data(select_bits,
			      multiple_of(l, BITS_PER_BYTE)/BITS_PER_BYTE);

	bitfield<uint64_t, uint64_t> walk_bits_bf;
	walk_bits_bf.set_raw_data(walk_bits, multiple_of(c*(l-1),
			   BITS_PER_TYPE(edge_datum_t))/BITS_PER_TYPE(edge_datum_t));

	// Infer the number of the starting vertex from the initial randomness
	if (sizeof(edge_datum_t) < sizeof(vertex_t)) {
		cerr << "(Lu bit extractor) Internal error: Assumption " <<
			"sizeof(edge_datum_t) < sizeof(vertex_t) failed!" << endl;
		exit(-1);
	}

	vertex_t vertex = *(vertex_t*)initial_rand;
	// Since w is aligned on sizeof(edge_datum_t), we can safely assume
	// that sizeof(vertex_t) bytes are available for the initial
	// vertex. We need, however, zero out the
	// bits [n_bits, sizeof(vertex_t)*BITS_PER_BYTE[ that are
	// not necessary for the required bit length.
	if (n_bits < sizeof(vertex_t)*BITS_PER_BYTE) {
		vertex &= ((vertex_t)1 << n_bits) - 1;
	}

	// After everything is set up, do the random walk.
	// In each step, the contribution to the parity result
	// is calculated.
	// TODO: Double-check all bit arithmetic operations in this part.
	edge_t next_edge;

	// Take the first of c walk steps, and compute l results this way
	uint64_t i, j;
	for (i = 0; i < c; i++) {
		// The vertex only changes the parity of the result if
		// the corresponding bit is set and select_bits_{i} == 1
		if (b.get_bit(vertex) && select_bits_bf.get_bit(i)) {
			res ^= 1;
//			cout << "    Flipping parity" << endl;
		}

		for (j = 0; j < l-1; j++) {
			// Choose the next edge from the initial randomness
			walk_bits_bf.get_bit_range((i*(l-1)+j)*bits_per_edge,
						   (i*(l-1)+j+1)*bits_per_edge-1,
						   reinterpret_cast<uint64_t*>(&next_edge));

			// ... and determine the number of the next vertex
			vertex = compute_next_vertex(vertex, next_edge);
		}
	}

	if (b.get_bit(vertex) && select_bits_bf.get_bit(i))
		res ^= 1;

	return res;
}
