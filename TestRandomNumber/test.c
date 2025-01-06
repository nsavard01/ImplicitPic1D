#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>
static uint64_t const multiplier = 6364136223846793005u;
static double const mult_factor = 1.0 / ((double)(UINT32_MAX) + 1.0);
static uint64_t const increment  = 1442695040888963407u;	// Or an arbitrary odd constant


double pcg32_random_r(int64_t* state)
{   
    uint64_t oldstate = (uint64_t)(*state);
    // Advance internal state
    *state = (int64_t)(oldstate * 6364136223846793005ULL + (increment|1));
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    uint32_t res = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    return (double)(res) * mult_factor;
}





