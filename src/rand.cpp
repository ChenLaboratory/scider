#include "rand.h"
#include "iostream"

//// PCG-PRNG by 2014 M.E. O'Neill / pcg-random.org
// This is the minimal pcg32 one-stream rng (64-bits states to 32-bits output) 
// with skip-ahead functionality.

// Seed the rng.
void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate)
{
  rng->state = 0U;
  pcg32_random_r(rng);
  rng->state += initstate;
  pcg32_random_r(rng);
}

// Generate a uniformly distributed 32-bit random number
uint32_t pcg32_random_r(pcg32_random_t* rng)
{
  uint64_t oldstate = rng->state;
  rng->state = oldstate * PCG_DEFAULT_MULTIPLIER_64 + PCG_DEFAULT_INCREMENT_64;
  uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
  uint32_t rot = oldstate >> 59u;
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// Skip-ahead the state of the RNG by delta in ln(delta) time complexity.
uint64_t pcg_advance_lcg_64(uint64_t state, uint64_t delta)
{
  uint64_t cur_mult = PCG_DEFAULT_MULTIPLIER_64;
  uint64_t cur_plus = PCG_DEFAULT_INCREMENT_64;
  uint64_t acc_mult = 1u;
  uint64_t acc_plus = 0u;
  while (delta > 0) {
    if (delta & 1) { // if delta is odd
      acc_mult *= cur_mult;
      acc_plus = acc_plus * cur_mult + cur_plus;
    }
    cur_plus = (cur_mult + 1) * cur_plus;
    cur_mult *= cur_mult;
    delta /= 2;
  }
  return acc_mult * state + acc_plus;
}


//// n choose k without replacement.
// This is a weird implementation that is highly specific to the task.
// Guarantees sampling k values within k picks (Important for stable 
// multi-threading since seed will be predictable).
// Samples are stored in the range x[n-k]...x[n-1] inclusive.
// x will be (partially) permuted after each run.
void sample_without_replacement(int* x,
                                int n,
                                int k,
                                pcg32_random_t* rng)
{
  int tempSwap;
  int j;
  double rng_val;
  for (int i=0;i<k;i++) {
    // Get a random 32-bits unsigned int then convert to range [0...1)
    // Conversion is slightly bias but fast & good enough when n << 2^32
    rng_val = (0x1.0p-32) * pcg32_random_r(rng); // [0...1)
    
    // Pick from n
    rng_val *= (--n); // [0...n-1)
    j = (int) rng_val;
    
    // Swap to put value in the correct place without removing stuffs from x 
    tempSwap = x[j];
    x[j] = x[n];
    x[n] = tempSwap;
  }
}