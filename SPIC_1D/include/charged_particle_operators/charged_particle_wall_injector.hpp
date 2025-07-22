
#pragma once
#include "charged_particle_operator.hpp"
#include "rand_gen/pcg_rng.hpp"
class charged_particle_wall_injector : public charged_particle_operator {

public:
    
    
    std::vector<double> current_density, v_x;
    std::vector<int> direction, wall_node_location;
    charged_particle_wall_injector(std::vector<int>& particle_indx, std::vector<double>& current_density, std::vector<double>& v_x, std::vector<int>& wall_node_location);
    void print_out() override;
    void run(const int thread_id, const double current_time, const double del_t, std::vector<charged_particle>& particle_list, const domain& world) override;

};


