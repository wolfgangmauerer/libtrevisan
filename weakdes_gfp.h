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

#ifndef WEAKDES_GFP_H
#define WEAKDES_GFP_H

#include "weakdes.h"
#include<vector>
#include<cmath>

class weakdes_gfp : public weakdes {
public:
	void compute_Si(uint64_t i, std::vector<uint64_t> &indices);
	uint64_t compute_d();
	void compute_admissible_params();

	long double get_r() { return 2*M_E ; }

private:
	unsigned int deg; // Degree of the polynomial used to compute (a, poly(a))

	template <class T, class F>
		T horner_poly(const std::vector<T> &coeffs, T x, T modulus, F inv_modulus);
};

#endif
