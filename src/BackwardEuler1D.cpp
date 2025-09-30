#include "BackwardEuler1D.h"
#include "TridiagonalSolver.h"

/**
 * @brief Construct a BackwardEuler1D solver.
 * 
 * Initializes the CFL number r and precomputes the tridiagonal system.
 *
 * The matrix A has the form:
 * 
 *   -r  1+2r  -r
 *   -r  1+2r  -r
 *        ...
 * 
 * which comes from discretizing (I - rA).
 */
BackwardEuler1D::BackwardEuler1D(
    const Grid1D& grid_,
    const HeatEquation1D& pde_,
    const BoundaryCondition& bc_,
    const InitialCondition& ic_, 
    double dt_, 
    double T_)
    : HeatSolver1D(grid_, pde_, bc_, ic_, dt_, T_)
{
    r = pde_.diffusion() * dt / (dx * dx);

    // The interior has Nx-1 unknowns excluding the boundaries. 
    u_interior.resize(Nx-1);

    // Precompute the tridiagonal matrix.
    lower.assign(Nx-2, -r);
    diag.assign(Nx-1, 1 + 2*r);
    upper.assign(Nx-2, -r);
}

/**
 * @brief Perform one Backward Euler time step.
 *
 * 1. Copy interior values from u into u_interior.
 * 2. Add contributions from boundary conditions to RHS.
 * 3. Solve the tridiagonal system (I - rA) u^{n+1} = u^n.
 * 4. Reapply Dirichlet boundary conditions.
 * 5. Copy updated interior solution back into u.
 *
 * @param t Current time
 */
void BackwardEuler1D::step(double t)
{
    // Step 1: Fill right-hand side with interior values of u.
    for (std::size_t i = 1; i < Nx; i++)
    {
        u_interior[i-1] = u[i];
    }

    // Step 2: Add contributions from Dirichlet boundary contributions
    u_interior[0]    += r * bc.left(t);   // Left boundary affects first row
    u_interior[Nx-2] += r * bc.right(t);  // Right boundary affects last row

    // Step 3: solve the tridiagonal system 
    tridiagonalSolve(lower, diag, upper, u_interior);

    // Step 4: Enforce Dirichlet boundary values on the entire vector.
    applyBoundaryConditions(t);

    // Step 5: Copy interior solutions back into u. 
    for (std::size_t i = 1; i < Nx; i ++)
    {
        u[i] = u_interior[i-1];
    }
}