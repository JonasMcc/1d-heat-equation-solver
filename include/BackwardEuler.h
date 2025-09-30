#pragma once

#include "HeatSolver1D.h"
#include <vector>

/**
 * @brief Backward Euler solver for the 1D heat equation.
 * 
 * Inherits from the HeatSolver1D and implements the implicit Backward Euler time step:
 * 
 *      (I - rA) u^{n+1} = u^n,
 * 
 * where r = α Δt / Δx² and A is the tridiagonal matrix representing the second derivative operator. 
 * 
 * @note Unlike Forward Euler, Backward Euler is unconditionally stable.
 */
class BackwardEuler1D : public HeatSolver1D 
{
public:
    /**
     * @brief Construct a BackwardEuler1D solver.
     * 
     * Precomputes the coefficients of the tridiagonal matrix equation.
     * 
     * @param grid Spatial grid
     * @param pde Heat equation parameters (diffusion)
     * @param bc Boundary condition
     * @param ic Initial condition
     * @param dt Time step size
     * @param T Final simulation time
     */
    BackwardEuler1D(const Grid1D& grid,
                    const HeatEquation1D& pde,
                    const BoundaryCondition& bc,
                    const InitialCondition& ic,
                    double dt,
                    double T);
protected:
    /**
     * @brief Perform a single Backward Euler time step. 
     * 
     * This function does the following:
     *    - Assembles the right-hand side of the equation.
     *    - Applies the boundary condition.
     *    - Solves the tridiagonal system
     *    - Updates the solution. 
     * 
     * @param t Current time
     */
    void step(double t) override;

private:
    double r;                                   //< CFL number
    std::vector<double> u_interior;             //< Solution on the interior.
    std::vector<double> lower, diag, upper;     //< The tridiagonal matrix
};