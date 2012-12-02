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

// A block weak design achieves r=1 by calling a primitive weak design in
// a special block mode

#include<iostream>
#include<cstdlib>
#include "weakdes_block.h"
#include "utils.hpp"
#ifndef NO_DEBUG
#include "debug.h"
#endif

using namespace std;

extern int debug_level;

void weakdes_block::compute_admissible_params() {
	// Since the underlying weak design might modify the parameters,
	// it must be registered before we can set the parameters.
	if (!wd) {
		cerr << "Internal error: Cannot set block design parameters "
		     << "without underlying weak design" << endl;
		exit(-1);
	}

	// t_requested and m_requested have been set by the call to set_params
	// that implicitely triggers this function
	wd->set_params(t_requested, m_requested);
	
	// The call to set_params has made sure that the parameters of
	// the underlying design have been adapted
	t = wd->get_t();
	m = wd->get_m();

	if (debug_level >= INFO) {
		cerr << "(Block design) Requested t: " << t_requested << endl;
		cerr << "(Block design) t compatible with constraints: " << t << endl;
	}
}

void weakdes_block::init_wd(block_t &blocks, class weakdes *wd) {
	this->blocks = blocks;
	this->wd = wd;
	initialised = true;
}

uint64_t weakdes_block::compute_d() {
	// The number of blocks l is given by blocks.nrow(),
	// and d=(l+1)*t^2
	if (!initialised) {
		cerr << "Internal error: Cannot compute d without given blocks" << endl;
		exit(-1);
	}

	return (blocks.nrow()+1)*t*t;
}

uint64_t weakdes_block::get_num_blocks() {
	if (!initialised) {
		cerr << "Internal error: Number of blocks requested before intialisation"
		     << endl;
		exit(-1);
	}

	return (blocks.nrow()+1);
}

void weakdes_block::compute_Si(uint64_t i, vector<uint64_t> &indices) {
	uint64_t local_i;
	tls_uint64_t::reference b_idx = tls_b_idx.local();
	tls_uint64_t::reference l_idx = tls_l_idx.local();
	tls_uint64_t::reference mult = tls_mult.local();
	tls_uint64_t::reference cached_i = tls_cached_i.local();
	tls_index_vec_t::reference cached_indices = tls_cached_indices.local();
	cached_indices.reserve(t);

	if (i == blocks(b_idx, HIGHEST_BIT)) {
		// We need to switch to the next block descriptor
		b_idx = b_idx + 1;
		l_idx = 0;
	}
		
	mult = l_idx % blocks(b_idx, BITS_PER_ITER);
	local_i = l_idx/blocks(b_idx, BITS_PER_ITER);

	if (blocks(b_idx, FIXUP) && mult == blocks(b_idx, BITS_PER_ITER)-1) {
		mult = blocks.nrow() - 1;
	}

	l_idx += 1;

	if (local_i == cached_i) {
		if (debug_level >= EXCESSIVE_INFO) {
			cerr << "Using cache for element " << local_i
			     << " (multiplication factor " << mult
			     << "*" << t*t << ")" << endl;
		}

		// Unsigned int is sufficient as data type because if the design
		// requires more initial randomness than that, it's pointless to
		// even try the computation.
		indices.clear();
		for (unsigned count = 0; count < t; count++) {
			indices.push_back(cached_indices[count] + mult*(t*t));
		}
	} else {
		// Call the basic construction
		wd->compute_Si(local_i, cached_indices);
		cached_i = local_i;

		copy_container(cached_indices, indices.begin());
	}
}
