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

#ifndef ONE_BITEXT_EXPANDER_H
#define ONE_BITEXT_EXPANDER_H

#include "1bitext.h"

typedef unsigned long long edge_datum_t;
typedef unsigned short edge_t;

class bitext_expander : public bitext {
public:
	bitext_expander(R_interp *r_interp) : bitext(r_interp) { };

	uint64_t get_l() const { return l; }
	uint64_t get_c() const { return c; }
	double get_w() const { return w; }

	void set_input_data(void *global_rand, struct phys_params& pp) override;

	// Pure virtual functions from the base class that need to be
	// implemented
	vertex_t num_random_bits();
	bool extract(void *rand);
	uint64_t compute_k() override;

private:
	void infer_params();

	// Helper functions
	vertex_t compute_next_vertex(vertex_t curr_vertex, edge_t next_edge);

	// TODO: These should probably go into a general helper class
	size_t multiple_of(size_t num, size_t mult);
	inline vertex_t_s do_mod(vertex_t_s a, vertex_t_s mod);

	double nu; // All parameters can be derived from w
	double w;
	vertex_t l, c; // c*(l-1) = number of random walk steps
	// The first of every c step is taken as result
	vertex_t sqrt_n; // n must be a power of two for this extractor
	unsigned short sqrt_n_bits; // How many bits are required to represent sqrt_n numbers
	unsigned short n_bits;      // Same for n

	// Mask to compute the upper and lower half of \mathbbm{Z}_{b}\otimes \mathbbm{Z}_{b}
	vertex_t zb_zb_mask;

	// Constants determined by the choice of underlying graph
	// NOTE: One of the assumptions of this code is that
	// g=2^{k} for some k.
	static const unsigned short g = 8; // g-regular graph
	static const unsigned short bits_per_edge = 3;
	static constexpr long double lambda0 = 5*sqrt(2.0)/8.0;

	// Offsets for the various portions of the initial randomness,
	// and the randomness overhead caused by alignment
	unsigned short additional_randomness ;
	unsigned short walk_bits_offset, select_bits_offset;
};

#endif
