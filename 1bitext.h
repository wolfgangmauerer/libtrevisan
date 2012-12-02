// -*- c++ -*-
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

#ifndef ONE_BITEXT_H
#define ONE_BITEXT_H

// NOTE: Unsigned long suffices for the number of vertex bits because 
// sizeof(unsigned long) is 32 resp 64 bits. Since very large numbers are
// supposed to be processed on 64 bit machines, we are restricted to
// 2^{62}/2^{3} = 2^{59} addressable bits, which is still considerably
// more than the amount of RAM contemporary machines can handle (2^{48}).
// However, we don't use mp_bitcnt_t as data type because the signed
// representation is required for the modulo calculations in the
// Gabber-Galil construction.
// NOTE: This definition would provide problems on Win64, where 
// unsigned long is 32 bits (we are safe on Unix/Linux).
typedef unsigned long vertex_t;
typedef signed long vertex_t_s;
#include<fstream>
#include<string>
#include<limits>
#include "bitfield.hpp"
#include "phys_params.h"
#include "R_interp.h"

class bitext {
public:
	bitext(R_interp *r_interp) {
		stat = nullptr;
		global_rand = nullptr;

		set_r_interp(r_interp);
	};

	virtual ~bitext() {
		if (stat)
			stat->close();
	};

	void set_stat_file(const std::string &filename) {
		// TODO: Open an output stream
	};

	virtual void set_input_data(void *global_rand, struct phys_params &pp) {
		this->pp = pp;
		this->global_rand = global_rand;
		mu = static_cast<long double>(pp.m)/(pp.alpha*pp.n);

		b.set_raw_data(global_rand, pp.n);
	};

	// Inform the bit extractor about the overlap parameter of the weak design
	//
	virtual void set_r(long double r) {
		this->r  = r;
	}

	// The following functions need to be implemented by derived
	// classes, that is, by specific 1-bit extractors.

	// Return the number of random input bits required for one
	// output bit
	virtual vertex_t num_random_bits() = 0;

	// Compute the required source entropy
	virtual uint64_t compute_k() = 0;

	// initial_rand: Pointer to randomness supplied by the weak design
	virtual bool extract(void *initial_rand) = 0;

private:
	void set_r_interp(R_interp *r_interp) {
		this->r_interp = r_interp;

		std::string bitext_code =
#include "generated/bitext_embedd.inc"
		r_interp->eval_code(bitext_code);
	};

	std::ofstream *stat;

protected:
	struct phys_params pp;
	long double mu;
	long double r;
	void *global_rand;
	bitfield<uint64_t, uint64_t> b;
	R_interp *r_interp;
};

#endif
