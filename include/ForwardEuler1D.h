#pragma once

#include "HeatSolver1D.h"
#include <vector>
#include <stdexcept>

/**
 * @brief Forward Euler solver for the 1D heat equation. 
 * 
 * Inherits from the HeatSolver1D and implements the Forward Euler time-step method: 
 * 
 *      u_i^{n+1} = u_i^n + r (u_{i-1}^n - 2 u_i^n + u_{i+1}^n)
 * 
 * where r = alpha * dt / dx^2. 
 * 
 * @note This method is conditionally stable. We require that r <= 0.5 for stability. 
 */
class ForwardEuler1D : public HeatSolver1D {
public:
    /**
     * @brief Construct a ForwardEuler1D solver.
     * 
     * @param grid Spatial grid
     * @param pde Heat equation parameters
     * @param bc Boundary conditions
     * @param ic Initial conditions. 
     * @param dt Time step size.
     * @param T Final simulation time. 
     * 
     * @throws std::invalid_argument if r > 0.5. 
     */
    ForwardEuler1D(const Grid1D& grid,
                   const HeatEquation1D& pde,
                   const BoundaryCondition& bc, 
                   const InitialCondition& ic,
                   double dt,
                   double T);

protected:
    /**
     * @brief Perform a single Forward Euler time step. 
     * 
     * Updates the solution vector u using the explicit finite difference formula.
     * 
     * @param t Current time. 
     */
    void step(double t) override; 

private: 
    double r;                           //< CFL number alpha*dt/dx^2
    std::vector<double> u_updated;      //< Temporary vector to store updated solution
};