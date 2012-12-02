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

#ifndef GF2ES_H
#define GF2ES_H
// Simple extension of NTL's GF2E class that allows for
// directly assigning unsigned long values without much ado

#include <NTL/GF2E.h>
#include <NTL/GF2X.h>
#include <NTL/GF2EX.h>
#include <NTL/tools.h>
#include <vector>

NTL_CLIENT

// "s" is for "sane value assignment"
// It's strange that this is seriously required -- seems
// fairly bogus. However, all conversion methods from long
// to GF2X and GF2E allow only for setting a _single_ bit,
// that is, the zeroeth coefficient. No idea what the rationale
// behind this is...
class GF2Es: public GF2E {
public:
    GF2Es() : GF2E() { }
    GF2Es(_ntl_ulong val) : GF2E() {
	(void)this->setValue(val);
    }

    GF2E *setValue(_ntl_ulong val) {
	if (val == 0) {
		clear(this->LoopHole());
		this->LoopHole().xrep[0] = 0;

		return this;
	}

	// Okay because _ntl_ulong is guaranteed to fit into one word
	this->LoopHole().xrep.SetLength(1);
	this->LoopHole().xrep[0] = val;
	rem(this->LoopHole(), this->LoopHole(), GF2E::modulus());

	return this;
    };

    GF2E *setValue(vector<_ntl_ulong> &val) {
	this->LoopHole().xrep.SetLength(val.size());

	bool all_zero = true;
	for (unsigned i = 0; i < val.size(); i++) {
		this->LoopHole().xrep[i] = val[i];
		if (val[i] != 0)
			all_zero = false;
	}

	if (all_zero) {
		clear(this->LoopHole());

		return this;
	}

	rem(this->LoopHole(), this->LoopHole(), GF2E::modulus());

	return this;
    };
};

// Same again for setting coefficients on GF2EX without having to explicitely
// create GF2X instances first.
inline void SetCoeff_s(GF2EX& x, long i, _ntl_ulong a) {
	if (a == 0)
		SetCoeff(x, i, GF2E::zero());
	else {
		GF2Es c(a);
		SetCoeff(x, i, c);
	}
}

#endif
