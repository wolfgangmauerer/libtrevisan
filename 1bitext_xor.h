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

#ifndef ONE_BITEXT_XOR_H
#define ONE_BITEXT_XOR_H

#include<limits>
#include "1bitext.h"
#include "bitfield.hpp"

// TODO: We could adaptivaly use a shorter datum if this has a performance benefit
typedef uint64_t datum_t;

class bitext_xor : public bitext {
public:
	bitext_xor(R_interp *r_interp) : bitext(r_interp) {
		l = std::numeric_limits<vertex_t>::max();
	};

	vertex_t get_l() const {
		return l;
	};

	void set_l(vertex_t l) {
		this->l = l;
	};

	void set_r(long double r) override {
		bitext::set_r(r);
		compute_l();
	};

	// Pure virtual functions from the base class that need to be
	// implemented
	vertex_t num_random_bits();
	bool extract(void *rand);
	uint64_t compute_k() override;

private:
	void compute_l();

	vertex_t l; // Number of random bits per run = l*log(n)
	bitfield<datum_t, uint64_t> local_bf;
};

#endif
