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

#ifndef WEAKDES_BLOCK_H
#define WEAKDES_BLOCK_H

#include<vector>
#include<limits>
#include<Rcpp.h>
#include<tbb/enumerable_thread_specific.h>
#include "weakdes.h"

typedef Rcpp::IntegerMatrix block_t;

// A block weak design is interesting from the object oriented design point of view:
// It IS-A weak design, and HAS-A weak design...
class weakdes_block : public weakdes {
public:
	weakdes_block() : tls_b_idx(static_cast<uint64_t>(0)),
			  tls_l_idx(static_cast<uint64_t>(0)),
			  // uint64_t_max means "cache invalid" by definition.
			  // This restricts the block design to a maximal bit length of
			  // 2^64-2 instead of 2^64-1. Poor us!
			  tls_cached_i(std::numeric_limits<uint64_t>::max()) {
		// NOTE: The argument-less superclass constructor is called
		// automatically
		wd = NULL;
		initialised = false;
	};

	long double get_r() { return 1.0 ; }
	uint64_t get_num_blocks();

	void init_wd(block_t &blocks, class weakdes *wd);
	void compute_Si(uint64_t i, std::vector<uint64_t> &indices);
	void compute_admissible_params();
	uint64_t compute_d();

private:
	bool initialised;

	block_t blocks;
	weakdes *wd;    // Basic construction the block design is based upon

	typedef tbb::enumerable_thread_specific<uint64_t> tls_uint64_t;
	typedef tbb::enumerable_thread_specific<std::vector<uint64_t>> tls_index_vec_t;
	tls_uint64_t tls_b_idx; // Block descriptor index
	tls_uint64_t tls_l_idx; // Local index: Bit within this descriptor range
	tls_uint64_t tls_mult; // Multiplier for set elements
	tls_uint64_t tls_cached_i;
	tls_index_vec_t tls_cached_indices;


// TODO: Turn into an enum
#define HIGHEST_BIT 0
#define ITER_PER_DESCRIPTOR 1
#define BITS_PER_ITER 2
#define FIXUP 3
};

#endif
