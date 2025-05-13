
#pragma once
#include <vector>
#include <string>
#include <memory>

class poisson_solver_1D {

protected:
    int number_unknowns; // number of unknowns in the system
    std::vector<double> source_term, solution; // diagonal and upper matrix elements
public:

    
    virtual ~poisson_solver_1D() = default;

    virtual void solve() = 0; // pure virtual function to solve the system

    inline const int& get_number_unknowns() const {
        return this->number_unknowns;
    };
    
    std::vector<double>& get_source_term() {
        return this->source_term;
    };

    std::vector<double>& get_solution() {
        return this->solution;
    };
};

