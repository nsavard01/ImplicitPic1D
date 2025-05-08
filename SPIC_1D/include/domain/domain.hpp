// include/domain.hpp
#pragma once
#include <vector>
#include <string>

class domain {

protected:
    size_t number_cells; // number of cells in the domain
    size_t number_nodes; // number of nodes in the domain
    double length_domain; // length of the domain
    double min_dx; // minimum cell size
    std::vector<double> grid_nodes; // grid nodes
    std::vector<double> cell_centers; // cell centers
    int right_boundary_condition; // right boundary condition type
    int left_boundary_condition; // left boundary condition type
public:

    
    virtual ~domain() = default;

    inline const double& get_min_dx() const {
        return this->min_dx;
    };

    inline const size_t& get_number_cells() const {
        return this->number_cells;
    };

    inline const size_t& get_number_nodes() const {
        return this->number_nodes;
    };

    inline const int& get_right_boundary_condition() const {
        return this->right_boundary_condition;
    };

    inline const int& get_left_boundary_condition() const {
        return this->left_boundary_condition;
    };

    // virtual void write_domain_to_file(const std::string& filename);
    // virtual void read_domain_from_file(const std::string& filename);
};


