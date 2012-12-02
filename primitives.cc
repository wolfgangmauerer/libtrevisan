#include<iostream>
#include "primitives.h"

using namespace std;

// Conversion routines between bit extractor/weak design constants and
// textual descriptions

wd_type get_weakdes(const string &wd) {
	if (wd == "gf2x")
		return wd_type::GF2X;

	if (wd == "gfp")
		return wd_type::GFP;

	if (wd == "block")
		return wd_type::BLOCK;

	if (wd == "aot")
		return wd_type::AOT;

	cerr << "Please specify a valid weak design construction!" << endl;
	exit(-1);
}

bext_type get_bitext(const string &bx) {
	if (bx == "xor")
		return bext_type::XOR;

	if (bx == "lu")
		return bext_type::LU;

	if (bx == "rsh")
		return bext_type::RSH;

	cerr << "Please specify a valid bit extractor construction!" << endl;
	exit(-1);
}

string bitext_to_string(bext_type bext) {
	switch (bext) {
	case bext_type::XOR:
		return "XOR";
		break;
	case bext_type::LU:
		return "Lu";
		break;
	case bext_type::RSH:
		return "Reed-Solomon-Hadamard";
		break;
	default:
		// Duh?!
		cerr << "Internal error: Unknown bit extractor specified" << endl;
		exit(-1);
	}
}

string weakdes_to_string(wd_type wdt) {
	switch (wdt) {
	case wd_type::GF2X:
		return "GF(2^m)";
		break;
	case wd_type::GFP:
		return "GF(p)";
		break;
	case wd_type::BLOCK:
		return "Block";
		break;
	case wd_type::AOT:
		return "Ahead of time";
		break;
	default:
		// Duh?!
		cerr << "Internal error: Unknown weak design specified" << endl;
		exit(-1);
	}
}
