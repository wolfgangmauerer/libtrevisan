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

// A weak design generator as described by Ma and Tan/Portmann
// Based on weakdes.cc. It remains to be seen if both variants
// should be merged back into a single implementation.
#include<sys/time.h>
#include<tr1/cstdint>
#include<cmath>
#include<string>
#include<iostream>
#include<cstdlib>
#include<bitset>
#include<vector>
#include<NTL/GF2E.h>
#include "GF2Es.h"
#include "irreps.h"
#include<NTL/GF2EX.h>
#include<NTL/GF2X.h>
#include<NTL/GF2XFactoring.h>
#include<NTL/tools.h>
#include "timing.h"
#include "utils.hpp"
#include "weakdes_gf2x.h"
#include "gf2n.h"
#ifndef NO_DEBUG
#include "debug.h"
#endif

#ifdef USE_NTL
NTL_CLIENT
#endif

using namespace std;

extern int debug_level;

#ifdef USE_NTL
void weakdes_gf2x::init_wd(unsigned int log_t) {
	// Set the minimal polynomial for the base extension field
	// GF(t) (where t = 2^k for some k)
	set_irrep(irred_poly, log_t);
	GF2E::init(irred_poly);
	
	if (debug_level >= INFO)
		cerr << "Picked weak design irreducible polynomial " << irred_poly << endl;
}
#else // !USE_NTL
void weakdes_gf2x::init_wd(unsigned int log_t) {
	// Choose the minimal polynomial from the precomputed list
	irred_poly = minweight_primpoly[log_t];

	// Perform the multiplications (after computing the helper mask)
	h = (1UL << (log_t-1));

	if (debug_level >= INFO)
		cerr << "Picked weak design irreducible polynomial " << irred_poly << endl;
}

// Horner's rule for polynomial evaluation, using GF(2^n) arithmetic
uint64_t weakdes_gf2x::horner_poly_gf2n(const vector<uint64_t> &coeffs, uint64_t x) {
	uint64_t res = 0;
 
	for(size_t i = 0; i < coeffs.size(); i++) {
	    res = gf2n_mult(res, x, irred_poly, h);
	    res = gf2n_add(res, coeffs[i]);
	}

	return res;
}
#endif

uint64_t weakdes_gf2x::compute_d() {
	return t*t;
}

void weakdes_gf2x::compute_admissible_params() {
	// This design has the constraint that t must be a power of two,
	// but works for arbitrary values of m

	// By default, assume that the parameters are alright.
	weakdes::compute_admissible_params();

	// Since the bit extractor requests t _bits_, it does not
	// matter if we work with a larger t, as long as we don't pass
	// more than t_requested indices per set upwards. This
	// requires more initial randomness than without constraints,
	// but does increase the number of loops to produce all output sets S_i.
	if (!(t_requested && !(t_requested & (t_requested - 1)))) {
		// Compute the closest power of 2 as effective t
		t = static_cast<uint64_t>(1) << ceil_log2(t_requested);

		cerr << "Warning (gf2x): t=" << t_requested << " is not in the form 2^k "
		     << "(using t=" << t << " instead)" << endl;
	}

	// Compute the degree of the evaluated polynomial
	// We need to use t_requested to determine the degree because with
	// a rounded t, we may end up in the case that m=t, which leads
	// to an impossible degree 0.
	deg = (unsigned int)ceil(log((long double)m)/log((long double)t_requested)-1);
	if (deg == 0) {
		cerr << "Error (gf2x): Polynomial degree is zero" << endl;
		cerr << "(m=" << m << ", t_requested=" << t_requested << ")" << endl;
		exit(-1);
	}

	if (debug_level >= INFO) {
		cerr << "gf2x-WD: log(m)=" << log((long double)m)
		     << ", log(t)=" << log((long double)t)
		     << ", polynomial degree=" << deg << endl;
	}
}

// We have a (t,r) weak design S_{1},...,S_{m} \in [d]
// with m is arbitrary, deg=ceil(log(m)/log(t)-1), d=t^{2},
// r=2e (e is Euler's constant), and t is a power of 2.

// This function computes the set S_{i}, which consists of
// t indices \in [d]. 
// NOTE: The effective parameters of the weak design
// are t and m, all others are dependent quantities.
// NOTE: In the end of the day, m denotes (after combining the weak
// design with an m-bit extractor) the number of output bits. See
// the comment in 1bitext.cc why uint64_t is more than sufficient 
// as data type for m, d, and t.
void weakdes_gf2x::compute_Si(uint64_t i, vector<uint64_t> &indices) {
	unsigned long count;
#ifdef USE_NTL
	GF2EX poly;
#else
	vector<uint64_t> coeff(deg+1);
#endif

#ifdef WEAKDES_TIMING
	meas_t start, end, delta;
	measure(&start);
#endif

#ifdef EXPENSIVE_SANITY_CHECKS
	if (!is_initialised()) {
		cerr << "Internal error: Weak design used before initalisation!" << endl;
		cerr << "(Did you call set_params()?)" << endl << flush;
		exit(-1);
	}
#endif

	// The field size is 2^l, which means we have values
	// from [0,2^m-1]. We therefore need to represent values
	// of at most t-1
	uint64_t log_t = numbits<uint64_t>(t-1);

	// Derive d and m from t and deg, ensuring that no
	// overflows can occur
	if (deg*numbits<uint64_t>(t) > sizeof(uint64_t)*8) {
	    cerr << "Internal error: Overflows when computing "
		 << "d=t^2 or m=t^deg" << endl;
	    exit(-1);
	}

	uint64_t d = t*t;
	
#ifdef EXPENSIVE_SANITY_CHECKS
	if (i > 0 && numbits<uint64_t>(i) > m) {
		cerr << "Internal error: i contains more than m bits!" << endl;
		exit(-1);
	}
#endif

	// i is interpreted as a log(m) bit quantity (since 0 <= i < m)
	// Compute the coefficients of the polynomial
	// by dividing the log(m) bits of i into deg parts with
	// log(t) bits each (deg = log(m)/log(t))

#ifdef EXPENSIVE_SANITY_CHECKS
	// We work in the extension field GF(2^t) and concatenate <a, poly(a)>.
	// Ensure that the concatenated result does not overflow the elementary
	// data type
	if (2*log_t > sizeof(uint64_t)*8) {
	    cerr << "Internal error: log_t is too large!" << endl;
	    exit(-1);
	}
#endif

	// Recall that t specifies the number of bit indices in S_{i},
	// m denote how many sets S_{i} are computed, and d
	// denotes the bit length of the input string (the random seed)
	if (debug_level >= EXCESSIVE_INFO) {
		cerr << "t=" << t << ", m=" << m << ", d=" << d
		     << ", log(t)=" << log_t << ", deg=" << deg << endl;
	}

	uint64_t mask = ((uint64_t)1 << log_t) - 1;

	// Set up a polynomial over GF(2^n) aka GF2[x], which is, in even
	// shorter terms, an element of GF(2^n)[x]
#ifdef USE_NTL
	uint64_t c;
	// NOTE: A polynomial of degree deg has deg+1 coefficients,
	// so we need to iterate up to and _including_ deg.
	for (count = 0; count <= deg; count++) {
		c = (i & (mask << (count*log_t))) >> count*log_t;

		if (debug_level >= EXCESSIVE_INFO)
			cerr << "Coefficient " << count << ": " <<
				c << endl;

		SetCoeff_s(poly, count, c);
	}

	if (debug_level >= INFO)
		cerr << "Constructed polynomial " << poly << endl;
#else
	// We need the coefficients for evaluating the polynomial
	// with Horner's rule, so pre-compute them in an extra loop
	for (count = 0; count < deg; count++) {
		// NOTE: Since the coefficients are restricted to
		// t bits, they are automatically members of GF(2^t)
		coeff.push_back((i & (mask << (count*log_t))) >> count*log_t);

		if (debug_level >= EXCESSIVE_INFO)
			cerr << "Coefficient " << count << ": " <<
				coeff[count] << endl;
	}
#endif

	// Iterate over all a \in GF(t) and compute <a, poly(a)>
	// NOTE: Since t fits into an elementary data type, poly(a)
	// is guaranteed to fit into an elementary data type when
	// evaluated over GF(t) as well (this need not be the case
	// if the polynomial were evaluated over \mathbb{R}).
	// NOTE: If the code turns out to be slow, optimisation would
	// most likely pay off here.
	uint64_t a;
	uint64_t res, res_num, res_num_tmp;

	// We only iterate up to t_requested and not up to the actual t
	// TODO: Ensure that this is not mathematically inadmissible for some reason.
	indices.clear();
	for (a = 0; a < t_requested; a++) {
#ifdef USE_NTL
		GF2Es val(a);
		GF2E res_gf2e;
		eval(res_gf2e, poly, val); // res = poly(val)
		if (IsZero(res_gf2e))
			res = 0;
		else
			res = res_gf2e.LoopHole().xrep[0];
#else
		res = horner_poly_gf2n(coeff, a);
#endif

		// Finally, generate the pair <a, poly(a)> as a concatenation
		// of the bit values, which denotes one element of S_{i}.
		// Since d=t^{2}, and we have made sure that d fits into
		// an elementary data type, the concatenation also fits into
		// an elementary data type.
#ifdef EXPENSIVE_SANITY_CHECKS
		size_t nb_a = numbits<uint64_t>(a);

		if (nb_a > log_t) {
			cerr << "Internal error: a contains more than " 
			     << "log(t)=" << log_t << " bits (" 
			     << nb_a << ")!" << endl;
			cout << "a = " << a << " bits: (" 
			     << bitset<64>(a) << ")" << endl;
			exit(-1);
		}
#endif

		res_num = a;
		res_num_tmp = res;

#ifdef EXPENSIVE_SANITY_CHECKS
		if ((res_num_tmp & mask) != res_num_tmp) {
		    cerr << "Internal error: res_num_tmp contains more "
			 << "than log_t bits" << endl;
		    exit(-1);
		}
#endif

		// Perform bitwise concatenation of a and poly(a)=res
		res_num_tmp <<= log_t;
		res_num |= res_num_tmp;

		// Although it is guaranteed that the concatenated
		// result fits into to available number of bits, the
		// numerical result may exceed d. This mandates a
		// modulo division (NOTE: Does not make any performance
		// difference if we use multiply_mod or the modulo
		// division here)
		// TODO: Not sure if this is correct for prime field
		// arithmetic, and double unsure if it's correct in GF(2^n).
		// Think of a generic method that allows to concatenate two
		// results without being prone to overflows.
		if (res_num > d)
		    res_num = res_num % d;

		// NOTE: Each result will be interpreted as a bit index
		// when the weak design is composed with a 1-bit extractor
		if (debug_level >= EXCESSIVE_INFO) {
		    cerr << "S contains <a, poly(a)>: <"
			 << bitset<15>(a) << ", " 
			 << bitset<15>(res) << ">" << endl;
			cerr << "Concatenated (i.e., bit index): "
			     << res_num << endl;
			cerr << "Binary: " << bitset<62>(res_num) << endl;
		}

		indices.push_back(res_num);
	} // End of iteration over all a \in GF(t_desired)

#ifdef WEAKDES_TIMING
	measure(&end);

	timestamp_subtract(&delta, &end, &start);
	if (stat != NULL)
	  (*stat) << delta_to_ns(delta) << endl;
#endif
}
