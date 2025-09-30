#pragma once 
#include <functional>

/**
 * @brief Represents Dirichlet boundary conditions for a 1D PDE. 
 * 
 * Stores the boundary values at the left and right ends of the domain.
 * The boundary conditions are functions of time: u(a, t) and u(b, t).  
 */
class BoundaryCondition {
public:
    /**
     * @brief Construct a BoundaryCondition object.
     * 
     * @param leftFunction Function defining the left boundary value u(a, t).
     * @param rightFunction Function defining the right boundary value u(b, t).  
     * 
     * @note The functions are expected to accept time (double) and return a double. 
     */
    BoundaryCondition(
        std::function<double(double)> leftFunction,
        std::function<double(double)> rightFunction);

    /**
     * @brief Evaluate the left boundary function at time t. 
     * 
     * @param t Time at which to evaluate the boundary function.
     * @return Value of left boundary function at time t. 
     */
    double left(double t) const; 

    /**
     * @brief Evaluate the right boundary function at time t. 
     * 
     * @param t Time at which to evaluate the boundary function.
     * @return Value of right boundary function at time t. 
     */
    double right(double t) const; 

private:
    std::function<double(double)> leftFunction;  //< Left boundary function
    std::function<double(double)> rightFunction; //< Right boundary function
};