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

#ifndef PRIMITIVES_H
#define PRIMITIVES_H

#include<string>
#include<cmath>
#include "1bitext_expander.h"
#include "1bitext_xor.h"
#include "1bitext_rsh.h"
#include "weakdes_gf2x.h"
#include "weakdes_gfp.h"
#include "weakdes_block.h"
#include "blockdes_params.h"
#include "weakdes_aot.h"

enum class wd_type { GF2X, GFP, BLOCK, AOT }; // Available weak designs
enum class bext_type { LU, XOR, RSH };   // Available one-bit extractors

// Map the weak design to the overlap parameter
// NOTE: While the aot weak design could have a different r than 2e, this
// does not make sense, so we ignore this case here -- it would lead to
// a causality dilemma in the main dispatcher, because, for instance, the
// XOR extractor needs r to compute l, and in turn the number of random bits
// per extraction run, while the weak designs depend on this information
// during their construction...
template<class T>
struct wd_overlap_trait {
	static constexpr long double r = 2*M_E;
};

template<>
struct wd_overlap_trait<class weakdes_block> {
	static constexpr long double r = 1.0;
};

// NOTE: We could also overload operator<< instead of providing string conversions
wd_type get_weakdes(const std::string &wd);
bext_type get_bitext(const std::string &bx);
std::string bitext_to_string(bext_type bext);
std::string weakdes_to_string(wd_type wdt);

#endif
