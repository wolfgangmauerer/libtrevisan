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

// Trevisan extractor main dispatcher
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring> // memset
#include <string>
#include <cmath>
#include <sys/time.h>
#include <tclap/CmdLine.h>
#include "ossl_locking.h"
#include "utils.hpp"
#include "debug_levels.h"
#include "timing.h"
#include "extractor.h"
#include "prng.hpp"
#include <tbb/tbb.h>
#include <tbb/task_scheduler_init.h>
#include <RInside.h>

using namespace TCLAP;
using namespace std;

int debug_level = 0;

// TODO: Use C++ memory management once the interferences between new and GMP
// are sorted out. Proper thing to do would be to only use the new operator,
// and write a failure handler (or use the provided exception)
void* do_alloc(size_t size, const string &type) {
	void *mem = new (nothrow) char[size];

	if (mem == NULL) {
		cerr << "Internal error: " <<
			"Cannot allocate space for " << type << endl;
		exit(-1);
	}

	return mem;
}

void *alloc_and_zero(size_t m) {
	int *mem = (int*)do_alloc((m/BITS_PER_TYPE(int)+1)*sizeof(int), "zero mem");

	memset(mem, 0, m/BITS_PER_BYTE+1);

	return mem;
}


// This function contains only generic calls to class bitext
// and weakdes. This makes the trevisan extractor algorithm
// independent of the specific 1-bit-extractor and weak design
// algorithms
inline bool trevisan_extract(uint64_t i, vector<uint64_t> &indices,
			     bitfield<unsigned int, uint64_t> &init_rand_bf,
			     unsigned int *y_S_i, params &params,
			     wd_file_lock_type &wd_file_lock, bitext *bext, weakdes *wd,
			     uint64_t t, int *extracted_rand, uint64_t num_rand_bits,
			     unsigned short offset) {
	if (debug_level >= PROGRESS) {
		if (i % 10000 == 0) {
			cout << "Computing bit " << i << endl;
		}
	}

	// Compute the index set S_{i}
	if (i % 10000 == 0 && debug_level >= INFO)
		cout << "Computing weak design for i=" << i << " (t="
		     << t << ")" << endl;

	wd->compute_Si(i, indices);

	if (params.save_weakdes) {
		{
			wd_file_lock_type::scoped_lock lock(wd_file_lock);
			params.wd_out.seekp(offset + i*t*sizeof(uint64_t), ios::beg);
			params.wd_out.write(reinterpret_cast<const char*>(&indices[0]),
					    t*sizeof(indices[0]));
		}
	}

	// Select the bits from y indexed by S_{i}, that
	// is, compute y_{S_{i}}. y_S_i thus contains t bits
	// of randomness.
	memset(y_S_i, 0, num_rand_bits/BITS_PER_BYTE+1);
	for (uint64_t count = 0; count < t; count++) {
		// Set the bit indexed count in y_S_i to the bit
		// indexed by indices[count] in init_rand.
#ifdef EXPENSIVE_SANITY_CHECKS
		if (indices[count] > num_rand_bits) {
			cerr << "Internal error: Bit index exceeds amount of "
			     << "available randomness (selected:" << indices[count]
			     << ", available: " << num_rand_bits
			     << ", index: " << count << ")!" << endl;
			exit(-1);
		}
#endif
		if (init_rand_bf.get_bit(indices[count])) {
			y_S_i[count/BITS_PER_TYPE(unsigned int)] |=
				((unsigned int)1 << (count % BITS_PER_TYPE(unsigned int)));
		}
	}

	// Feed the subset of the initial randomness to the
	// 1-bit extractor, and insert the resulting bit into
	// position i of the extracted randomness.
	if (!params.skip_bitext) {
		if(bext->extract(y_S_i)) {
			extracted_rand[i/BITS_PER_TYPE(int)] |=
				(1 << (i % BITS_PER_TYPE(int)));
		}
	}
}

void trevisan_dispatcher(class weakdes *wd, class bitext *bext, params &params) {
	// To determine the amount of required randomness, first determine
	// how many bits the single bit extractor requires per extracted bit.
	// Then, pass this value to the weak design, which infers
	// the total amount d of requires random input bits.
	// NOTE: For the 1-bit-extractor run, the number of required
	// random bits is given by num_rand_bits() -- it
	// can be slightly larger than the theoretically required
	// amount because of alignment or efficiency reasons.
	uint64_t t = bext->num_random_bits();
	cout << "1-bit extractor requires " << t << " bits per "
	     << "extracted bit" << endl;
	wd->set_params(t, params.pp.m);

	uint64_t num_rand_bits = wd->compute_d();
	long double r = wd->get_r();
	if (params.verbose) {
		cout.precision(3);
		cout << "Total amount of initial randomness: " << num_rand_bits
		     << " (" << fixed << (((float)num_rand_bits)/(1024*1024)) << " MiBit)"
		     << endl;
		cout << "Ratio of extracted bits to initial randomness: "
		     << ((long double)params.pp.m)/num_rand_bits << endl;
		cout << "-------------------------------------------------------------"
		     << endl << endl;
	}

	if (params.dryrun)
		return;

	unsigned short offset = 0;
	if (params.save_weakdes) {
		// Save information about the weak design characteristics
		// in a "header"
		// The first entry specifies the starting offset for
		// the actual weak design sets
		offset = sizeof(unsigned short) +
			sizeof(wd_type) + 3*sizeof(uint64_t);
		params.wd_out.write((char*)&offset, sizeof(unsigned short));
		params.wd_out.write((char*)&(params.wdt), sizeof(wd_type));
		params.wd_out.write((char*)&(params.pp.m), sizeof(uint64_t));
		params.wd_out.write((char*)&t, sizeof(uint64_t));
		params.wd_out.write((char*)&num_rand_bits, sizeof(uint64_t));
		params.wd_out.write((char*)&r, sizeof(long double));
	}

	// Allocate space for the extracted randomness
	int *extracted_rand = (int*)alloc_and_zero(params.pp.m);
	
	// Create the initial randomness y (and make the result accessible
	// as void pointer)
	vector<uint64_t> initial_rand_vector;
	// TODO: Using uint64_t instead of the previously used, but incorrect
	// unsigned int as index data type induces a 20% performance penalty.
	// We should thus adaptively use shorter data types for all bit index related
	// operations when smaller amounts of data are processed.
	bitfield<unsigned int, uint64_t> init_rand_bf;
	init_rand_bf.set_raw_data(create_randomness(num_rand_bits, initial_rand_vector),
				  num_rand_bits);

	wd_file_lock_type wd_file_lock;

	// TODO: Adapt this to the generic bitext/weakdesign API
	ofstream *bitext_stat_file;
	ofstream *wd_stat_file;
	// TODO: Let the user select this at run time
	// TODO: Access to the statistics files must be serialised
	if (0) {
	  bitext_stat_file = new ofstream();
	  bitext_stat_file->open("/tmp/bitext_stats.txt");

	  wd_stat_file = new ofstream();
	  wd_stat_file->open("/tmp/wd_stats.txt");
	} else {
	  bitext_stat_file = NULL;
	  wd_stat_file = NULL;
	}

	meas_t start, end, delta;
	measure(&start);

	// 3.) Compute all weak design index sets S_{i}, select the
	// appropriate bits from the initial randomness, and pass the
	// selected randomness to the 1-bit extractor. Concatenate
	// the 1-bit results.
	parallel_for(tbb::blocked_range<uint64_t>(0, params.pp.m),
		     [bext, wd, t, extracted_rand, num_rand_bits, offset,
		      &init_rand_bf, &params, &wd_file_lock]
		     (const tbb::blocked_range<uint64_t>& range) {
			     vector<uint64_t> indices;
			     vector<unsigned int> y_S_i;

			     y_S_i.reserve(num_rand_bits/BITS_PER_TYPE(unsigned int)+1);
			     indices.reserve(t);

			     // TODO: Create an apply class and move all constant
			     // parameters to private data elements that are initialised
			     // in the constructor. Essentially, the only parameter
			     // that should remain in the function call is i
			     for (auto i = range.begin(); i != range.end(); i++) {
				     trevisan_extract(i, indices, init_rand_bf, &y_S_i[0],
						      params, wd_file_lock, bext, wd, t,
						      extracted_rand, num_rand_bits, offset);
			     }
		     });

	measure(&end);

	timestamp_subtract(&delta, &end, &start);

	cout.precision(3);
	cout << endl << "-------------------- Summary --------------------" << endl;
	cout << "Required " << delta_to_ms(delta) << " ms for "
	     << params.pp.m << " bits extracted from " << params.pp.n << " bits." << endl;

	cout << "Performance: " << params.pp.m/(delta_to_s(delta)*1000)
	     << " kbits/s" << endl;

	cout << "Ratio of extracted bits to initial randomness: "
	     << ((long double)params.pp.m)/num_rand_bits << endl;


	if (bitext_stat_file) {
	  bitext_stat_file->close();
	  wd_stat_file->close();
	}

	// TODO: Save the extracted randomness to a file
}

void parse_cmdline(struct params &params, int argc, char **argv) {
	try {
		CmdLine cmd("Trevisan extractor framework", ' ', "Sep 2012");

		ValueArg<unsigned long> nArg("n", "inputsize",
					"Length of input data (bits)",
					true, 0, "Integer");
		ValueArg<unsigned long> mArg("m", "outputsize",
					"Length of extracted data (bits)",
					true, 0, "Integer");
		ValueArg<string> weakdesArg("w", "weakdes",
					    "Weak design construction "
					    "(gf2x, gfp, block, aot)",
					   false, "gf2x", "string");
		ValueArg<string> basic_weakdesArg("", "basic-weakdes",
						  "Basic weak design construction for the "
						  "block weak design (gf2x, gfp, aot)",
						  false, "gf2x", "string");
		ValueArg<string> bitextArg("x", "bitext",
					    "Bit extractor construction (lu, xor, rsh)",
					    false, "xor", "string");
		ValueArg<double> alphaArg("a", "alpha",
					  "Source entropy factor alpha (0 < alpha < 1)",
					  false, 0.9, "Real");
		ValueArg<double> epsArg("e", "eps", "1 bit error probability (0 < eps < 1)",
					false, 1e-7, "Real");
		ValueArg<double> lu_nuArg("", "lu:NU",
					       "nu parameter for the Lu extractor",
						false, 0.45, "Integer");
		ValueArg<int> numTasksArg("", "numtasks",
					  "Number of tasks (0 means unlimited)",
					  false, -1, "Integer");
		SwitchArg ignoreEntropyArg("", "ignore-entropy-violation",
					   "Allow minimum entropy requirements",
					   false);
		SwitchArg bytesArg ("", "bytes", "Use bytes instead of bits "
				    "to compute input/output sizes", false);
		SwitchArg kiloArg ("", "kilo", "Multiply size units by 1024", false);
		SwitchArg megaArg ("", "mega", "Multiply size units by 1024*1024", false);
		SwitchArg verboseArg ("v", "verbose", "Enable verbose output", false);
		SwitchArg dryrunArg ("", "dry-run", "Only compute parameters, don't extract",
				     false);

		cmd.add(nArg);
		cmd.add(mArg);
		cmd.add(weakdesArg);
		cmd.add(basic_weakdesArg);
		cmd.add(alphaArg);
		cmd.add(epsArg);
		cmd.add(lu_nuArg);
		cmd.add(numTasksArg);
		cmd.add(bitextArg);
		cmd.add(ignoreEntropyArg);
		cmd.add(bytesArg);
		cmd.add(kiloArg);
		cmd.add(megaArg);
		cmd.add(verboseArg);
		cmd.add(dryrunArg);

		SwitchArg skipBitextArg ("", "skip-bitext", "Skip the bit extraction step",
					 false);
		cmd.add(skipBitextArg);

		ValueArg<string> weakdesFileArg("", "weakdes-file",
						"File to save/load the weak design",
						false, "", "string");
		cmd.add(weakdesFileArg);

#ifndef NO_DEBUG
		MultiSwitchArg debugArg ("d","debug",
				      "Emit diagnostic output "
				      "(use multiple times for more details)");
		cmd.add(debugArg);
#endif

		cmd.parse(argc, argv);

		params.pp.n = nArg.getValue();
		params.pp.m = mArg.getValue();

		params.wdt = get_weakdes(weakdesArg.getValue());
		if (params.wdt == wd_type::BLOCK) {
			params.basic_wdt = get_weakdes(basic_weakdesArg.getValue());
			if (params.basic_wdt == wd_type::BLOCK) {
				cerr << "Cannot use the block design as basic design for "
				     << "block design" << endl;
				exit(-1);
			}
		}

		params.bxt = get_bitext(bitextArg.getValue());

		unsigned int mult_factor = 1;
		if (kiloArg.getValue() && megaArg.getValue()) {
			cout << "Please specify either kilo or mega as multiplier, "
			     << "not both!" << endl;
			exit(-1);
		}

		if (bytesArg.getValue()) {
			mult_factor *= 8;
		}
		if (kiloArg.getValue()) {
			mult_factor *= 1024;
		}
		if (megaArg.getValue()) {
			mult_factor *= 1024*1024;
		}

		params.pp.n *= mult_factor;
		params.pp.m *= mult_factor;

		params.verbose = verboseArg.getValue();

		// Dry run does not make any sense without verbose
		// parameter messages
		params.dryrun = dryrunArg.getValue();
		if (params.dryrun)
			params.verbose = true;

		params.ignore_entropy = ignoreEntropyArg.getValue();

		params.pp.alpha = alphaArg.getValue();
		if (params.pp.alpha <= 0 || params.pp.alpha >= 1) {
			cerr << "Source entropy factor alpha must be in the range (0,1)"
			     << endl;
			exit(-1);
		}

		params.pp.eps = epsArg.getValue();
		if (params.pp.eps <= 0 || params.pp.eps > 1) {
			cerr << "1 bit error eps must be in the range (0,1)"
			     << endl;
			exit(-1);
		}

		params.pp.lu_nu = lu_nuArg.getValue();
		if (params.pp.lu_nu <= 0 || params.pp.lu_nu >= 0.5) {
			cerr << "nu parameter for the Lu extractor must satisfy "
			     << "0 < nu < 0.5" << endl;
			exit(-1);
		}

		params.num_tasks = numTasksArg.getValue();
#ifndef NO_DEBUG
		debug_level = debugArg.getValue();
#endif

		params.skip_bitext = skipBitextArg.getValue();
		params.wd_filename = weakdesFileArg.getValue();
		params.save_weakdes = false;

		if (params.wdt == wd_type::AOT ||
		    (params.wdt == wd_type::BLOCK && params.basic_wdt == wd_type::AOT)) {
			// When the AOT weak design was specified, the filename
			// is not used to save, but to read the data
			if (params.wd_filename == "") {
				cerr << "Cannot use AOT weak design without filename "
				     << "(specify --weakdes=file=<file>)" << endl;
				exit(-1);
			}
		} else if (params.wd_filename != "") {
			// No AOT design used, but weak design filename given
			// -> Save the weak design
			params.save_weakdes = true;

			params.wd_out.open (params.wd_filename, ios::out | ios::binary);
			if(!params.wd_out.is_open()) {
				cout << "Error: Cannot open weak design output file!"
				     << endl;
				exit(-1);
			}
		}
	} catch (ArgException &e) {
		cerr << "Error encountered during command line parsing: "
		     << e.error() << " for argument " << e.argId() << endl;
		exit(-1);
	}
}


// Some weak designs need special intialisation sequences
// This can be handled by providing an appropriate template
// specialisation for the initialisation function that does nothing
// in the general case
template<class W>
void init_primitives(W *wd, class bitext *bext, struct params &params) {
	bext->set_r(wd_overlap_trait<W>::r);
	return;
}

template<>
void init_primitives(class weakdes_gf2x *wd, class bitext *bext, struct params &params) {
	// Initialise the weak design GF(2^t)
	// There are 2^t field elements, which means we have values
	// from [0,2^t-1]. We therefore need to represent values
	// of at most t-1 bits.
	bext->set_r(wd_overlap_trait<weakdes_gf2x>::r);
	uint64_t t = bext->num_random_bits();
	int log_t = numbits<uint64_t>(t-1);
	wd->init_wd(log_t);
}

template<>
void init_primitives(class weakdes_aot *wd, class bitext *bext, struct params &params) {
	bext->set_r(wd_overlap_trait<weakdes_aot>::r);
	wd->set_file(params.wd_filename);
	wd->set_file_lock(new wd_file_lock_type);
}


// When there's a sub-design, the main design can only be a block design
template<class S>
void init_primitives(class weakdes_block *wd, S *sub_wd, class bitext *bext,
		     struct params &params) {
	blockdes_params bd_params(params.R);
	block_t blocks;

	bext->set_r(wd_overlap_trait<weakdes_block>::r);
	uint64_t t = bext->num_random_bits();
	blocks = bd_params.compute_blocks(2*M_E, params.pp.m, t);
	wd->init_wd(blocks, sub_wd);

	return;
}

void alloc_init_weakdes(wd_type type, weakdes **wd, bitext *bext, params &params) {
	switch (type) {
	case wd_type::GF2X: {
		*wd = new weakdes_gf2x;
		init_primitives(dynamic_cast<weakdes_gf2x*>(*wd), bext, params);
		break;
	}
	case wd_type::GFP: {
		*wd = new weakdes_gfp;
		init_primitives(dynamic_cast<weakdes_gfp*>(*wd), bext, params);
		break;
	}
	case wd_type::AOT: {
		*wd = new weakdes_aot;
		init_primitives(dynamic_cast<weakdes_aot*>(*wd), bext, params);
		break;
	}
	default:
		cerr << "Internal error: Unknown weak design type requested" << endl;
		exit(-1);
	}
}


void show_params(struct params &params, bitext *bext, weakdes *wd) {
	cout << "-------------------------------------------------------------"
	     << endl;
	cout << "Number of input bits: " << params.pp.n << endl;
	cout << "Number of extracted bits: " << params.pp.m << endl;
	cout << "Bit extractor: " << bitext_to_string(params.bxt) << endl;
	cout << "Weak design: " << weakdes_to_string(params.wdt);
	if (params.wdt == wd_type::BLOCK) {
		cout << " (basic construction: "
		     << weakdes_to_string(params.basic_wdt) << ", number of blocks: "
		     << dynamic_cast<weakdes_block*>(wd)->get_num_blocks()
		     << ")" << endl;
	}
	cout << endl;

	cout << "Source entropy factor alpha: " << params.pp.alpha << endl;
	cout << "1 bit error eps: " << params.pp.eps << endl;
	cout << "Required source entropy: " << bext->compute_k() << " (available: "
	     <<  static_cast<uint64_t>(params.pp.alpha*params.pp.n) << ")" << endl;
	if (params.bxt == bext_type::XOR) {
		cout << "XOR-extractor parameter l: "
		     << dynamic_cast<bitext_xor*>(bext)->get_l() << endl;
	}

	if (params.bxt == bext_type::LU) {
		cout << "Lu-extractor parameters: nu=" << params.pp.lu_nu
		     << ", c=" << dynamic_cast<bitext_expander*>(bext)->get_c()
		     << ", l=" << dynamic_cast<bitext_expander*>(bext)->get_l()
		     << ", w=" << dynamic_cast<bitext_expander*>(bext)->get_w() << endl;
	}
	if (params.save_weakdes && params.verbose)
		cout << "Saving weak design to " << params.wd_filename << endl;
	cout << "Using " << params.num_tasks << " parallel computation unit";
	if (params.num_tasks != 1)
		cout << "s";
	cout << endl;
}


/////////////////////////////////// Dispatcher /////////////////////////////
int dispatch(struct params &params) {
	int num_tasks;
	num_tasks = tbb::task_scheduler_init::default_num_threads();
	if (params.num_tasks > 0)
		num_tasks = params.num_tasks;
	else
		params.num_tasks = num_tasks;
	tbb::task_scheduler_init init(num_tasks);

	// We rely on openssl being compiled with thread support
#define OPENSSL_THREAD_DEFINES
#include <openssl/opensslconf.h>
#ifndef OPENSSL_THREADS
#error Please use a thread capable openssl
#endif

	init_timekeeping();
	init_openssl_locking();
	// NOTE: Must be set up before any Rcpp data types are instantiated
	R_interp *r_interp = new R_interp;

	// Suitable quick test parameters for XOR and GF2X: eps=1e-7, n=10e9, m=10e5
	// Parameters of Ma et al.: RSH and GF2X, n=2**14, m=2**13
	class bitext *bext;
	class weakdes *wd;

	// We rely on double having at least 64 bits in the modular
	// multiplication code
	if (sizeof(double) < 8) {
	    cerr << "Internal error: double must at least encompass 64 bits!"
		 << endl;
	    exit(-1);
	}

	vector<uint64_t> global_rand_vector;
	void *global_rand = create_randomness(params.pp.n, global_rand_vector);

	switch (params.bxt) {
	case bext_type::LU: {
		bitext_expander *btx_expander = new bitext_expander(r_interp);
		btx_expander->set_input_data(global_rand, params.pp);
		bext = btx_expander;
		break;
	}
	case bext_type::XOR: {
		bitext_xor *btx_xor = new bitext_xor(r_interp);
		btx_xor->set_input_data(global_rand, params.pp);
		bext = btx_xor;
		break;
	}
	case bext_type::RSH: {
		bitext_rsh *btx_rsh = new bitext_rsh(r_interp);
		btx_rsh->set_input_data(global_rand, params.pp);
		bext = btx_rsh;
		break;
	}
	default:
		cerr << "Internal error: Unknown 1-bit extractor requested" << endl;
		exit(-1);
	}

	if (params.wdt == wd_type::BLOCK) {
		weakdes *sub_wd;
		alloc_init_weakdes(params.basic_wdt, &sub_wd, bext, params);

		wd = new weakdes_block;
		init_primitives(dynamic_cast<weakdes_block*>(wd), sub_wd, bext, params);
	} else {
		alloc_init_weakdes(params.wdt, &wd, bext, params);
	}

	if (params.verbose)
		show_params(params, bext, wd);

	// Make sure that the source provides a sufficient amount of entropy
	if (bext->compute_k() > params.pp.alpha*params.pp.n) {
		if (params.ignore_entropy)
			cerr << "Warning: ";
		else
			cerr << "Error: ";
		cerr << "Source does not contain sufficient entropy "
		     << "for specified extraction parameters!" << endl;
		if (!params.ignore_entropy) {
			cerr << "(Choose a smaller m, change eps or use a "
			     << "different weak design to satisfy the constraints," << endl
			     << "or specify --ignore-entropy-violation to "
			     << "ignore the constraint)" << endl;
			exit(-1);
		}
	}

	// Finally, call the trevisan extractor proper
	trevisan_dispatcher(wd, bext, params);

	if (params.save_weakdes)
		params.wd_out.close();

	delete bext;
	delete wd;
	delete r_interp;

	return 0;
}

int main(int argc, char** argv) {
	struct params params;
	parse_cmdline(params, argc, argv);
	int ret = dispatch(params);
	
	return(ret);
}
