// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
// Slight modification where it passes reference to int64_t for interoperability with fortran, and returns a double precision float

#include <stdint.h>
static double const mult_factor_PCG = 1.0 / ((double)(UINT32_MAX) + 2.0); // Multiplier for conversion to [0,1] double


double pcg32_random_r(int64_t* state)
{   
    // Convert state to uint64_t
    uint64_t oldstate = (uint64_t)(*state);
    // Advance internal state
    *state = (int64_t)(oldstate * 6364136223846793005ULL + 1442695040888963407ULL);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    uint32_t res = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    return ((double)(res) + 1.0) * mult_factor_PCG;
}





