#ifndef BITFIELD_H
#define BITFIELD_H

// A bitfield based on raw data considerably larger than the largest elementary type
// with the ability to select sub-ranges of bits
// This might seem like a good use case for a STL bitfield, but the differences
// are large enough to make a specific implementation worth while.

#include<tr1/cstdint>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cstring> // memset
#include<bitset>
#include "debug_levels.h"

extern int debug_level;

// C is the type used for chunking; I is the index type
template<typename C, typename I>
class bitfield {
public:
	bitfield() : data(nullptr), n(0) {
		bits_per_type=sizeof(*data)*8;
	};

	bitfield(void *data, I n) {
		bits_per_type=sizeof(*data)*8;
		set_raw_data(data, n);
	};

	bitfield(C *data, I n) {
		bits_per_type=sizeof(*data)*8;
		this->data = data;
		this->n = n;
	};

	void set_raw_data(void *data, I n) {
		this->data = reinterpret_cast<C*>(data);
		this->n = n;
	};

	// Caller is responsible to provide enough storage space in res
	// Select bit range [start, end]
	// NOTE: Bit indices are, as usual, zero-based
	void get_bit_range(I start, I end, C *res) {
#ifdef EXPENSIVE_SANITY_CHECKS
		if (end > n || start > n || start >= end) {
			// TODO: Throw an exception instead
			std::cerr << "Internal error: Bit index bogous" << std::endl;
			std::cerr << "(start=" << start << ", end=" << end
				  << ", n=" << n << ")" << std::endl;
			std::exit(-1);
		}
#endif

		memset(res, 0, ceil(((double)(end - start + 1))/8));
		// Compute chain elements and bit within this element for the
		// start and end positions
		I start_chunk = start/bits_per_type;
		I end_chunk = end/bits_per_type;

		unsigned short start_bit = start % bits_per_type;
		unsigned short end_bit;

		I dest_chunk = 0;
		C chunk;
		C mask;

		if (start_chunk == end_chunk)
			end_bit = end % bits_per_type;
		else
			end_bit = bits_per_type - 1;

		mask = compute_mask(start_bit, end_bit);
		res[dest_chunk] = (data[start_chunk] & mask) >> start_bit;

		if (start_chunk == end_chunk)
			return;

		I dest_start_bit = end_bit - start_bit + 1;
		if (dest_start_bit == bits_per_type) {
			dest_chunk++;
			dest_start_bit = 0;
		}

		I dest_bits;

		for (I curr_chunk = start_chunk+1; curr_chunk < end_chunk; curr_chunk++) {
			// For the inner chunks, we can always select the full chunk from
			// the input data, but need to split it across the output
			// data field
			chunk = data[curr_chunk];

			// How many bits remain in the destination chunk
			dest_bits = bits_per_type - dest_start_bit;

			// Fill up the current destination chunk
			mask = compute_mask(0, dest_bits - 1);
			res[dest_chunk] |= ((chunk & mask) << dest_start_bit);
			dest_chunk++;

			// ... and fill the next destination chunk as far as
			// possible unless the previous chunk was completely
			// drained and there is nothing left for the new chunk
			if (dest_bits != bits_per_type) {
				mask = compute_mask(dest_bits, bits_per_type - 1);
				res[dest_chunk] = (chunk & mask) >> dest_bits;
			}

			// Compute new starting position in the destination chunk
			dest_start_bit = bits_per_type - dest_bits;
		}

		end_bit = end % bits_per_type;
		dest_bits =  bits_per_type - dest_start_bit;
		
		mask = compute_mask(0, end_bit);
		chunk = data[end_chunk] & mask;

#ifdef EXPENSIVE_SANITY_CHECKS
		if (debug_level >= EXCESSIVE_INFO) {
			std::cerr << "end_bit: " << end_bit << ", dest_bits: " << dest_bits
				  << ", mask: " << std::bitset<32>(mask) << ", end_chunk: "
				  << end_chunk << std::endl;
			std::cerr << "data[end_chunk]: " << std::bitset<32>(data[end_chunk])
				  << std::endl;
			std::cerr << "masked data:     " << std::bitset<32>(chunk)
				  << std::endl;
			std::cerr << "dest_start_bit: " << dest_start_bit << std::endl;
		}
#endif
		// Any excess bits that do not fit into the current result chunk
		// are shifted out to the left here
		res[dest_chunk] |= (chunk << dest_start_bit);

		if (end_bit + 1 > dest_bits) {
			// We need to split the final chunk contribution
			// across two result chunks
			dest_chunk++;

			// Shift out the already consumed bits to the right,
			// and include the remaining bits into the final result chunk
			res[dest_chunk] = (chunk >> dest_bits);
		} 
	};

	inline bool get_bit(I i) {
#ifdef EXPENSIVE_SANITY_CHECKS
		if (i > n) {
			// TODO: Throw an exception instead
			std::cerr << "Internal error: Bit index out of range" << std::endl;
			std::cerr << "(total length: " << n << ", requested: "
				  << i << ")" << std::endl;
			std::exit(-1);
		}
#endif

		// Compute chain element the bit is contained in
		I idx = i/bits_per_type;
		C chunk = data[idx];

		bool bit = chunk & (static_cast<I>(1) << (i % bits_per_type));

		return bit;
	};

private:
	C *data;
	uint64_t n; // Length of the field in bits
	unsigned short bits_per_type;

	inline C compute_mask(I start_bit, I end_bit) {
		C mask = 0;
#ifdef EXPENSIVE_SANITY_CHECKS
		assert(end_bit < bits_per_type);
		assert(end_bit >= start_bit);
#endif
		if (end_bit < bits_per_type-1)  {
			mask = (1 << (end_bit+1)) - 1;
		} else {
			mask = ~static_cast<C>(0);
		}

		mask &= ~((static_cast<C>(1) << start_bit) - static_cast<C>(1));
		
		return mask;
	}
};

#endif
