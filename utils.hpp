#ifndef UTILS_HPP
#define UTILS_HPP

// don't use namespace std -- we may want to use this in code
// where std is not imported into the main namespace

#include<cstddef>
#include<cassert>
#include<cstdlib>
#include<cmath> // floating point logarithm
#include<iostream>
#include<limits>

#define BITS_PER_BYTE	8
#define BITS_PER_TYPE(_type)	(sizeof(_type)*BITS_PER_BYTE)

// Compute \lceil log_{2}(n)\rceil. We don't care about efficiency here.
template<class T>
std::size_t ceil_log2(T num) {
	assert(num > 0);

	if (num & ((T)1 << (BITS_PER_TYPE(T) - 1))) {
		if ((num & ((T)1 << (BITS_PER_TYPE(T) - 1) - 1)) > 0) {
			std::cerr << "Internal error: Overflow in floor_log2" 
				  << std::endl;
			std::exit(-1);
		}

		return BITS_PER_TYPE(T);
	}

	std::size_t log;
	for (log = BITS_PER_TYPE(T)-1; log--; log >= 0) {
		if (num & ((T)1 << log))
			if ((num & (((T)1 << log) - 1)) > 0)
				return log+1;
			else
				return log;
	}

	// Should never be reached
	assert(1==0);
	return 42;
}

// Determine how many bits are required to store a number num
// Could as well compute ceil_log2<T>(num+1)
template<class T>
std::size_t numbits(T num) {
	assert(num >= 0);
	std::size_t bits;

	if (num == 0)
	    return 1;

	bits = ceil_log2<T>(num);

	if (num & ((T)1 << bits)) {
		return bits + 1;
	}

	return bits;
}

// Same for floor.
template<class T>
std::size_t floor_log2(T num) {
	assert(num > 0);
	size_t log;

	for (log = BITS_PER_TYPE(T)-1; log--; log >= 0) {
		if (num & ((T)1 << log))
			return log;
	}
}

// Shannon entropy
template<class T>
T h(T x) {
	if (x < 2*std::numeric_limits<T>::epsilon())
		return(0);

	if (x > 1-2*std::numeric_limits<T>::epsilon())
		return(0);

	return(-x*log2(x) -(1-x)*log2(1-x));
}

template<class C, class OutIter>
OutIter copy_container(const C& c, OutIter result) {
	return std::copy(c.begin(), c.end(), result);
}

#endif
