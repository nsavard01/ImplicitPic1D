// include/domain.hpp
#pragma once
#include <vector>
#include <string>
#include <memory>

class domain {
    
public:

    int number_cells; // number of cells in the domain
    int domain_type;
    int number_nodes; // number of nodes in the domain
    double length_domain; // length of the domain
    double min_dx; // minimum cell size
    std::vector<double> dx_dxi; // cell sizes in the domain
    std::vector<double> grid_nodes; // grid nodes
    std::vector<double> cell_centers; // cell centers
    int right_boundary_condition; // right boundary condition type
    int left_boundary_condition; // left boundary condition type
    virtual ~domain() = default;

    inline const int& get_domain_type() const {
        return this->domain_type;
    };

    inline const double& get_min_dx() const {
        return this->min_dx;
    };

    inline const std::vector<double>& get_dx_dxi() const {
        return this->dx_dxi; // return the cell sizes
    }

    inline const int& get_number_cells() const {
        return this->number_cells;
    };

    inline const double& get_domain_length() const {
        return this->length_domain;
    };

    inline const int& get_number_nodes() const {
        return this->number_nodes;
    };

    inline const std::vector<double>& get_grid() const {
        return this->grid_nodes; // return the grid nodes
    };

    inline const int& get_right_boundary_condition() const {
        return this->right_boundary_condition;
    };

    inline const int& get_left_boundary_condition() const {
        return this->left_boundary_condition;
    };


    virtual void print_out() = 0;
    // virtual void write_domain_to_file(const std::string& filename);
    // virtual void read_domain_from_file(const std::string& filename);
};

std::unique_ptr<domain> create_domain_from_file(const std::string& filename, int scheme_type);
