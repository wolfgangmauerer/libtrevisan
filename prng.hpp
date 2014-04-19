#ifndef PRNG_H
#define PRNG_H

#include <random>
#include <vector>
#include <cstdint>
#include <functional>

template<class T>
void *create_randomness(uint64_t n, std::vector<T> &rand) {
	// Prepare n bits of randomness using a pseudo-random number generator.
	// A vector instead of a raw array is used to store the data because
	// this way, it's easier for the caller to keep track of the allocation.
	std::uniform_int_distribution<T>
		distribution(0, std::numeric_limits<T>::max());
	std::mt19937 mt_engine(42);
	auto generator = std::bind(distribution, mt_engine);
	
	uint64_t num_chunks = n/sizeof(rand[0])+1;
	rand.reserve(num_chunks);
	
	for (auto i = 0; i < num_chunks; i++) {
		rand[i] = generator(); 
	}

	return (&rand[0]);
}

#endif
