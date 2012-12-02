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

#ifndef WEAKDES_H
#define WEAKDES_H

#include<fstream>
#include<string>
#include<vector>
#include<tr1/cstdint>
#include <tbb/spin_mutex.h>

// Lock type is typedef'd so that it can be exchanged in a single spot as
// per the TBB recommendation
typedef tbb::spin_mutex wd_file_lock_type;

class weakdes {
public:
	weakdes() {
#ifdef EXPENSIVE_SANITY_CHECKS
		initialised = false;
#endif
		stat = NULL;
	};   
	
	weakdes(const std::string filename) {
		set_stat_file(filename);
	};
	
	virtual ~weakdes() {
		if (stat)
			stat->close();
	};
	
	void set_stat_file(const std::string &filename) {
		// TODO: Open an output stream
	};
	
	// NOTE: Checking for an initialised design is expensive because we
	// have to do it in compute_Si, so the check happens in the inner loop
	bool is_initialised() {
#ifdef EXPENSIVE_SANITY_CHECKS
		return initialised;
#else
		return true;
#endif
	};

	
	virtual void set_params(uint64_t t, uint64_t m) {
		t_requested = t;
		m_requested = m;

		compute_admissible_params();
#ifdef EXPENSIVE_SANITY_CHECKS
		initialised = true;
#endif
	};	

	virtual void compute_admissible_params() {
		// By default, we assume that the weak design implementation
		// works for every desired parameter combination. Variants
		// for which this is not the case need to overwrite this
		// function.
		t = t_requested;
		m = m_requested;
	}
	
	// NOTE: The caller is responsible that indices contains
	// enough space for t instances of uint64_t
	virtual void compute_Si(uint64_t i, std::vector<uint64_t> &indices) = 0;

	virtual uint64_t compute_d() = 0;
	virtual long double get_r() = 0; // Overlap factor

	uint64_t get_t() { return t; };
	uint64_t get_m() { return m; };
	
protected:
	uint64_t t; // Elements per output set, that is, |S_i|
	uint64_t m; // Number of output sets, that is, |{S_{i}}|

	// Depending on the implementation of the weak design, there can be
	// constraints on the admissible values of t. Consequently, we store two
	// sets of values: The desired values of t and degree, and the admissible
	// values of t and degree compute_admissible_params() can be provided by
	// weak design implementations to fill in the admissible values given the
	// desired values.
	uint64_t t_requested;
	uint64_t m_requested;

	std::ofstream *stat;
	
private:
#ifdef EXPENSIVE_SANITY_CHECKS
	bool initialised;
#endif
};

#endif
