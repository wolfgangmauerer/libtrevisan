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

#ifndef R_INTERP_H
#define R_INTERP_H

#include<RInside.h>
#include<memory>
#include<tbb/spin_mutex.h>

typedef tbb::spin_mutex r_interp_lock_t;

class R_interp {
public:
	R_interp() {
		if (++instances == 1) {
			r_interp.reset(new RInside(0, NULL));
		}
	};

	void eval_code(const std::string &code);
	int parse_eval(const std::string &call, SEXP &ans);

private:
	// Make sure there are no concurrent calls to the R sesson
	static r_interp_lock_t r_lock;
	static std::unique_ptr<RInside> r_interp;
	static short instances; // Ensure there's only one R session

	// Ensure that the object cannot be copied or assigned
	R_interp& operator = (const R_interp& other) { }
	R_interp(const R_interp& other) { }
};

#endif
