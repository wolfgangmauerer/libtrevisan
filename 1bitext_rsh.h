// -*- C++ -*-
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

#ifndef ONE_BITEXT_RSH_H
#define ONE_BITEXT_RSH_H

#include "1bitext.h"
#include "bitfield.hpp"
#include<vector>
#ifdef USE_NTL
#include<NTL/GF2EX.h>
#else
#include<openssl/bn.h>
#endif

#ifdef USE_NTL
NTL_CLIENT
#endif

class bitext_rsh : public bitext {
public:
	bitext_rsh(R_interp *r_interp) : bitext(r_interp) { };
	~bitext_rsh();

	void set_input_data(void *global_rand, struct phys_params &pp) override {
		bitext::set_input_data(global_rand, pp);
		compute_r_l();

		// The coefficients of the polynomial (computed from the
		// global randomness) do not vary between invocations, so we
		// can compute them before the extraction starts.
		create_coefficients();
	};

	// Pure virtual functions from the base class that need to be
	// implemented
	vertex_t num_random_bits();
	bool extract(void *initial_rand);
	uint64_t compute_k() override;
	
private:
	void create_coefficients();
	void compute_r_l();

	uint64_t r;      // Order of polynomial, see the paper for details
	uint64_t l;      // Seed length parameter, see the paper for details
	uint64_t chars_per_half; // Number of char instances in half of the intial randomness

	typedef uint64_t idx_t; // Index type for the bit field

#ifdef USE_NTL
	GF2EX poly;
	GF2X irred_poly;
	typedef _ntl_ulong chunk_t; // Chunk type for the bitfield
	typedef _ntl_ulong data_t;  // Internal library type to store GF2m elements
#else
	std::vector<BIGNUM*> coeffs;
	int irred_poly[5]; // The irrep is either a trinomial or pentanomial
	BIGNUM *horner_poly_gf2n(BIGNUM *x);

	// TODO: Determine if smaller index types cause any significant speedup
	typedef uint64_t chunk_t;
	typedef BN_ULONG data_t;
#endif

	bitfield<chunk_t, idx_t> b;
};

#endif
