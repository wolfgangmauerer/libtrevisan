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

#include "weakdes_aot.h"
#include "primitives.h"
#include <iostream>
#ifndef NO_DEBUG
#include "debug.h"
#endif

using namespace std;

extern int debug_level;

void weakdes_aot::set_file(const string &filename) {
	wd_type wd;
	dat.open (filename.c_str(), ios::in | ios::binary);
	if(!dat.is_open()) {
		cout << "Error: Cannot open weak design input file!"
		     << endl;
		exit(-1);
	}

	dat.seekg (0, ios::beg);
	dat.read((char*)&offset, sizeof(unsigned short));
	dat.read((char*)&wd, sizeof(wd_type));
	dat.read((char*)&m, sizeof(uint64_t));
	dat.read((char*)&t, sizeof(uint64_t));
	dat.read((char*)&d, sizeof(uint64_t));
	dat.read((char*)&r, sizeof(long double));

	if (debug_level >= PROGRESS) {
		cerr << "Reading precomputed weak design from " << filename
		     << " with parameters t=" << t << ", m=" << m << ", d=" << d
		     << ", r=" << r << ", data offset=" << offset << endl;
		cerr << "Original weak design type: " << weakdes_to_string(wd) << endl;
	}
}

uint64_t weakdes_aot::compute_d() {
	return d;
}

long double weakdes_aot::get_r() {
	if (!dat.is_open()) {
		cerr << "(aot weak design) Internal error: Cannot determine r without "
		     << "input file!" << endl;
		exit(-1);
	}

	return r;
}

void weakdes_aot::compute_Si(uint64_t i, vector<uint64_t> &indices) {
	// TODO: This is a hot-spot, use a lock-less construction in
	// the read path.
	if (file_lock != nullptr)
		flock.acquire(*file_lock);

	dat.seekg(offset + i*t*sizeof(uint64_t), ios::beg);
	// Since the array is continuous, we can read all indices in one single go
	dat.read((char*)indices[0], sizeof(uint64_t)*t);

	if (file_lock != nullptr)
		flock.release();
}

void weakdes_aot::compute_admissible_params() {
	weakdes::compute_admissible_params();

	// The provided value for t can be larger than the requested
	// value, see the comment in weakdes_gf2x.cc. But not the other way.
	if (t_requested > t) {
		cerr << "Error (weakdes_aot): Requested t=" << t_requested << " is larger "
		     << "than the stored value ot t (" << t << "), aborting." << endl;
		exit(-1);
	}
}
