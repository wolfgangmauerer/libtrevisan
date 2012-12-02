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

#include "tbb/spin_mutex.h"
#include "R_interp.h"

using namespace std;

short R_interp::instances = 0;
r_interp_lock_t R_interp::r_lock;
std::unique_ptr<RInside> R_interp::r_interp;

void R_interp::eval_code(const string &code) {
	{
		r_interp_lock_t::scoped_lock(r_lock);
		r_interp->parseEval(code);
	}
}

int R_interp::parse_eval(const string &call, SEXP &ans) {
	int res;
	{
		r_interp_lock_t::scoped_lock(r_lock);
		res = r_interp->parseEval(call, ans);
	}

	return res;
}
