#include "ForwardEuler1D.h"

/**
 * @brief Construct a ForwardEuler1D solver.
 *
 * Initializes the solver by computing the CFL number r = alpha*dt/dx^2
 * and allocating a temporary vector for the updated solution.
 *
 * @param grid_ Spatial grid
 * @param pde_ Heat equation parameters (diffusion coefficient α)
 * @param bc_ Boundary conditions
 * @param ic_ Initial condition
 * @param dt_ Time step size
 * @param T_  Final simulation time
 *
 * @throws std::invalid_argument if the time step is too large (r > 0.5),
 *         which would make the Forward Euler method unstable.
 */
ForwardEuler1D::ForwardEuler1D(
    const Grid1D& grid_,
    const HeatEquation1D& pde_,
    const BoundaryCondition& bc_,
    const InitialCondition& ic_,
    double dt_, 
    double T_)
    : HeatSolver1D(grid_, pde_, bc_, ic_, dt_, T_)
{
    // Compute CFL number
    r = pde.diffusion() * dt / (dx * dx);

    // Allocate temporary vector to store updated solution. 
    u_updated.resize(grid.size(), 0.0);

    // Stability check for forward Euler:
    if (r > 0.5)
    {
        throw std::invalid_argument("Unstable time step: r > 0.5.");
    }
}

/**
 * @brief Perform one Forward Euler time step.
 *
 * Updates the solution vector u using the explicit finite difference formula:
 *      
 *      u_i^{n+1} = u_i^n + r * (u_{i-1}^n - 2*u_i^n + u_{i+1}^n)
 * 
 * Boundary conditions are applied to both the current and updated solution vectors.
 *
 * @param t Current simulation time
 */
void ForwardEuler1D::step(double t)
{
    // Apply boundary conditions to the current solution
    applyBoundaryConditions(t);

    // Update interior points using Forward Euler finite difference method. 
    for(std::size_t i = 1; i < Nx; i++)
    {
        u_updated[i] = u[i] + r * (u[i-1] - 2*u[i] + u[i+1]);
    }

    // // Apply boundary conditions to the updated vector
    applyBoundaryConditionsUpdated(t, u_updated);

    // Copy updated solution back to main vector
    u = u_updated;
}

