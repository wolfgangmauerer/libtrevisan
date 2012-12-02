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

// A weak design generator as described by Hartman and Raz
// Not (yet) for distribution
#include<sys/time.h>
#include<tr1/cstdint>
#include<cmath>
#include<string>
#include<iostream>
#include<cstdlib>
#include<bitset>
#include<vector>
#include<gmp.h>
#include "timing.h"
#include "utils.hpp"
#include "weakdes_gfp.h"
#ifndef NO_DEBUG
#include "debug.h"
#endif

using namespace std;

extern int debug_level;

// Algorithm taken from Matters Computational, 39.1.2.
// Works for moduli of up to 61 bits.
// NOTE: We need a 64 bit float type, so double (instead
// of long double as used by fxt) suffices on 64 bit machines
// NOTE: It's _extremely_ important to inline this function -- gives
// a 35% performance increase on my machine.
inline uint64_t multiply_mod(uint64_t a, uint64_t b, uint64_t m,
			     double m_inv) {
    uint64_t x = a*b;
    uint64_t y = m*(uint64_t)((double)a*(double)b*m_inv +
			      (double)1/2);
    uint64_t r = x-y;

    if ((int64_t)r < 0)
	r += m;

    return r;
}

// Horner's rule for polynomial evaluation, using modular arithmetic
template <class T, class F>
T weakdes_gfp::horner_poly(const vector<T> &coeffs, T x, T modulus, F inv_modulus) {
	T res = 0;
 
	for(size_t i = 0; i < coeffs.size(); i++) {
	    res = multiply_mod(res, x, modulus, inv_modulus) + coeffs[i];
	}

	// Since the coefficients are upper bounded by the modulus,
	// the final non-modular addition can make the result too large
	// by at most modulus-2. Consequently, subtracting the modulus
	// once to get the remainder in case the result is too large
	// suffices (brings about 10% performance improvement
	// over multiply_mod(res, 1, modulus, inv_modulus))
	if (res >= modulus)
	    res -= modulus;

#ifdef EXPENSIVE_SANITY_CHECKS
	if (res >= modulus) {
	    cerr << "Internal error: assumption in horner_poly failed!"
		 << endl;
	    exit(-1);
	}
#endif

	return res;
}

void weakdes_gfp::compute_admissible_params() {
	// This design has the constraint that t must be a prime _number_
	// (it could also work with arbitrary prime powers, but we don't
	// support this -- use the GF(2^n) specialisation in this case),
	// but works for arbitrary values of m

	// By default, assume that the parameters are alright.
	weakdes::compute_admissible_params();

	// The primality test used in the following is probabilistic.
	// With more iterations, the results become more certain
	static unsigned int test_repetitions = 10;

	// If not, fix up
	mpz_t t_gmp;
	mpz_init_set_ui(t_gmp, t_requested);

	if (mpz_probab_prime_p(t_gmp, test_repetitions) != 2) {
		// Pick next prime that is larger than t_requested
		mpz_t new_prime; mpz_init(new_prime);
		mpz_nextprime(new_prime, t_gmp);

		if (mpz_fits_ulong_p (new_prime) != 0) {
			t = mpz_get_ui(new_prime);
		} else {
			// Okay, it's essentially impossible that the
			// new prime does not fit into an unsigned
			// long, but let's do our due diligence
			cerr << "Internal error: New prime does not fit into elementary "
			     << " data types!" << endl;
			exit(-1);
		}

		cerr << "Warning: t is not a prime number "
		     << "(using t=" << t << " instead)" << endl;
	}

	// Compute the degree of the evaluated polynomial
	// See weakdes_gf2x for why we use t_requested.
	deg = (unsigned int)ceil(log((long double)m)/log((long double)t_requested)-1);

	if (deg == 0) {
		cerr << "Error (gfp): Polynomial degree is zero" << endl;
		cerr << "(m=" << m << ", t_requested=" << t_requested << ")" << endl;
		exit(-1);
	}

	if (debug_level >= INFO) {
		cerr << "gfp-WD: log(m)=" << log((long double)m)
		     << ", log(t)=" << log((long double)t)
		     << ", polynomial degree=" << deg << endl;
	}
}

uint64_t weakdes_gfp::compute_d() {
	return t*t;
}

// TODO: Update the code to account for the fact that we can have
// arbitrary values of m that are not subject to any condition
// TODO: This should also simplify the parameter estimation.

// We have a (t,r) weak design S_{1},...,S_{m} \in [d]
// with m is arbitrary, t=p^{deg}, (deg \in \mathbbm{N}), d=t^{2},
// r=2e (e is Euler's constant), and p prime.

// This function computes the set S_{i}, which consists of
// t indices \in [d]. 
// NOTE: The effective parameters of the weak design
// are t and deg, the others are dependent quantities.
// (naturally, other combinations can also be used
// as defining parameters)
// NOTE: In the end of the day, m denotes (after combining the weak
// design with an m-bit extractor) the number of output bits. See
// the comment in 1bitext.cc why uint64_t is more than sufficient 
// as data type for m, d, and t.
void weakdes_gfp::compute_Si(uint64_t i, vector<uint64_t> &indices) {
	unsigned long count;
	vector<uint64_t> coeff(deg+1);

#ifdef WEAKDES_TIMING
	meas_t start, end, delta;
	measure(&start);
#endif

#ifdef EXPENSIVE_SANITY_CHECKS
	// Derive d and m from t and deg, ensuring that no
	// overflows can occur
	if (deg*numbits<uint64_t>(t) > sizeof(uint64_t)*8) {
	    cerr << "Internal error: Overflows when computing "
		 << "d=t^2 or m=t^deg" << endl;
	    exit(-1);
	}
#endif

	uint64_t d = t*t;
	uint64_t m = (uint64_t)pow((long double)t, (long double)deg);
	
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
	uint64_t log_t = numbits<uint64_t>(t);

	// We work mod t and concatenate <a, poly(a)>. Ensure that the
	// concatenated result does not overflow the elementary data type

#ifdef EXPENSIVE_SANITY_CHECKS
	if (2*log_t > sizeof(uint64_t)*8) {
	    cerr << "Internal error: log_t is too large!" << endl;
	    exit(-1);
	}
#endif

	// Recall that t specifies the number of bit indices in S_{i},
	// m denote how many sets S_{i} are computed, and d
	// denotes the bit length of the input string (the random seed)
	if (debug_level >= EXCESSIVE_INFO) {
		cerr << "t=" << t << ", m=" << m << ", d=" << d << endl;
		cerr << "log(t)=" << log_t << endl;
	}

	uint64_t mask = ((uint64_t)1 << log_t) - 1;

	// We need the coefficients for evaluating the polynomial
	// with Horner's rule, so pre-compute them in an extra loop
	// NOTE: A polynomial of degree deg has deg+1 coefficients,
	// so we need to iterate up to and _including_ deg.
	for (count = 0; count <= deg; count++) {
		coeff[count] = (i & (mask << (count*log_t))) >> count*log_t;
		coeff[count] = coeff[count] % t;

		if (debug_level >= EXCESSIVE_INFO)
			cerr << "Coefficient " << count << ": " <<
				coeff[count] << endl;
	}

	// Iterate over all a \in GF(t) and compute <a, poly(a)>
	// NOTE: Since t fits into an elementary data type, poly(a)
	// is guaranteed to fit into an elementary data type when
	// evaluated over GF(t) as well (this need not be the case
	// if the polynomial were evaluated over \mathbb{R}).
	// NOTE: If the code turns out to be slow, optimisation would
	// most likely pay off here.
	uint64_t a;
	uint64_t res, res_num, res_num_tmp;
	double inv_t; // Precomputed floating-point 1/t
	
	// We only iterate up to t_requested and not up to the actual t
	// TODO: Ensure that this is not mathematically inadmissible for some reason.
	indices.clear();
	for (a = 0; a < t_requested; a++) {
		inv_t = (double)1/t;
		res = horner_poly<uint64_t, double>(coeff, a, t, inv_t);

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
		if (res_num > d)
		    res_num = res_num % d;

		// NOTE: Each result will be interpreted as a bit index
		// when the weak design is composed with a 1-bit extractor
		if (debug_level >= EXCESSIVE_INFO) {
		    cerr << "S contains <a, poly(a)>: <"
			 << bitset<30>(a) << ", " 
			 << bitset<30>(res) << ">" << endl;
			cerr << "Concatenated (i.e., bit index): "
			     << res_num << endl;
			cerr << "Binary: " << bitset<62>(res_num) << endl;
		}

		indices.push_back(res_num);
	}

#ifdef WEAKDES_TIMING
	measure(&end);

	timestamp_subtract(&delta, &end, &start);
	if (stat != NULL)
	  (*stat) << delta_to_ns(delta) << endl;
#endif
}
