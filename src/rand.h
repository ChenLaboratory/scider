#ifndef RAND_H
#define RAND_H

#include <cstdint> //for uint64_t
#include <Rinternals.h>
#include <R_ext/Random.h>

// PCG-PRNG
#define PCG_DEFAULT_MULTIPLIER_64  6364136223846793005ULL
#define PCG_DEFAULT_INCREMENT_64   1442695040888963407ULL

struct pcg_state_64 {    // Internals are *Private*.
  uint64_t state;             // RNG state.  All values are possible.
};
typedef struct pcg_state_64 pcg32_random_t;

void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate);
uint32_t pcg32_random_r(pcg32_random_t* rng);
uint64_t pcg_advance_lcg_64(uint64_t state, uint64_t delta);

// Sampling
void sample_without_replacement(int* x, int n, int k, pcg32_random_t* rng);

#endif