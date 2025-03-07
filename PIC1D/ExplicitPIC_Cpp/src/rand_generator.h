#pragma once

#include <cstdint>
#include <omp.h>
#include <ctime>

static const double mult_factor_PCG = 1.0 / (static_cast<double>(UINT32_MAX) + 2.0); // Multiplier for conversion to [0,1] double


    
double pcg32_random_r(uint64_t* state)
{
    // Convert state to uint64_t
    uint64_t oldstate = *state;
    // Advance internal state
    *state = oldstate * 6364136223846793005ULL + 1442695040888963407ULL;
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = static_cast<uint32_t>((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = static_cast<uint32_t>(oldstate >> 59u);
    uint32_t res = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    return ((static_cast<double>(res)) + 1.0) * mult_factor_PCG;
}
