#include <iostream>
#include <stdio.h>
#include "pcg_rng.h"
#include <omp.h>
#include <iomanip>
#include <cmath>
#include "Constants.h"
#include "global_inputs.h"
#include "domain.h"
#include "particle.h"


int main() {
    int thread_id, i;
    int number_stats = 10;
    double sum = 0.0;
    double var = 0.0;
    double val;
    number_nodes = 101;
    
    Domain world("../InputData/Geometry.inp");
    read_particle_inputs("../InputData/ParticleTypes.inp");
    omp_set_num_threads(2);
    #pragma omp parallel private(thread_id, i, val) reduction(+:sum, var)
    {   
        thread_id = omp_get_thread_num();
        // // Initialize PRNG state for each thread
        initialize_pcg(false); 
        // cg_state = static_cast<uint64_t>(12345 + thread_id);
        // #pragma omp critical
        // {
           
        // }
        
        
        #pragma omp critical
        {
            for (int i = 0; i < number_stats; i++) {
            std::cout << "Thread " << thread_id << " has val " << pcg32_random_r() << std::endl;
    
        }
        
        }
        
    }

    return 0;
}
