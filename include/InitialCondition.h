#pragma once
#include <functional>

/**
 * @brief Represents the initial condition of a 1D PDE. 
 * 
 * Stores a function of space x that defines the initial condition u(x, 0). 
 */
class InitialCondition {
public: 
    /**
     * @brief Constrct a InitialCondition object. 
     * 
     * @param f Function defining the initial condition u(x, 0).
     * 
     * @note The function is expected to accept a spatial coordinate x (double) and return a double. 
     */
    InitialCondition(std::function<double(double)> f);

    /**
     * @brief Evaluate the initial condition at a given x. 
     * 
     * @param x The given spatial coordinate. 
     * @return The value of the initial condition evaluated at x.
     */
    double operator()(double x) const; 
private:
    std::function<double(double)> f; //< Function storing the initial condition. 
};