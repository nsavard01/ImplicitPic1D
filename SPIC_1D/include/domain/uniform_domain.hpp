// include/domain.hpp
#pragma once
#include "domain/domain.hpp"
#include <vector>
#include <iostream>

class uniform_domain : public domain {

public:

    
    uniform_domain(size_t num_cells, double length_domain, int left_boundary_condition, int right_boundary_condition) {
        this->number_cells = num_cells;
        this->length_domain = length_domain;
        this->left_boundary_condition = left_boundary_condition;
        this->right_boundary_condition = right_boundary_condition;
        this->min_dx = length_domain / num_cells;
        this->number_nodes = num_cells + 1;
        this->grid_nodes.resize(this->number_nodes);
        this->cell_centers.resize(this->number_cells);
        for (size_t i = 0; i < this->number_nodes; i++) {
            this->grid_nodes[i] = i * this->min_dx;
        }
        for (size_t i = 0; i < this->number_cells; i++) {
            this->cell_centers[i] = (this->grid_nodes[i] + this->grid_nodes[i + 1]) * 0.5;
            printf("Cell center %zu: %f\n", i, this->cell_centers[i]);
        }
    }


};
