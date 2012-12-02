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

#ifndef WEAKDES_AOT_H
#define WEAKDES_AOT_H

// A weak design computed ahead of time
#include<vector>
#include<string>
#include<fstream>
#include "weakdes.h"

class weakdes_aot : public weakdes {
public:
	weakdes_aot() : file_lock(nullptr) {};

	~weakdes_aot() {
		if (dat.is_open())
			dat.close();
	};

	void set_file(const std::string &filename);
	void set_file_lock(wd_file_lock_type *file_lock) {
		// Synchronises file output
		this->file_lock = file_lock;
	};
	void compute_Si(uint64_t i, std::vector<uint64_t> &indices);
	uint64_t compute_d();
	void compute_admissible_params();
	long double get_r();

private:
	std::ifstream dat;
	uint64_t d;
	long double r;
	unsigned short offset;
	wd_file_lock_type *file_lock;
	wd_file_lock_type::scoped_lock flock;
};

#endif
