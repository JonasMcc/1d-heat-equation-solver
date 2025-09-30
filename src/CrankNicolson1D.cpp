#include "CrankNicolson1D.h"
#include "TridiagonalSolver.h"

/**
 * @brief Construct a CrankNicolson1D solver.
 *  
 * Initializes the CFL number r and sets up the constant tridiagonal system (I - r/2 A). 
 *
 * LHS matrix stencil:
 *   -r/2   1+r   -r/2
 */
CrankNicolson1D::CrankNicolson1D(
    const Grid1D& grid_,
    const HeatEquation1D& pde_,
    const BoundaryCondition& bc_,
    const InitialCondition& ic_, 
    double dt_,
    double T_)
    : HeatSolver1D(grid_, pde_, bc_, ic_, dt_, T_)
{
    r = pde_.diffusion() * dt / (dx * dx);

    // Precompute constant LHS system. 
    lower.assign(Nx - 2, -r / 2.0);
    diag.assign(Nx - 1, 1.0 + r);
    upper.assign(Nx - 2, -r / 2.0);

    rhs.resize(Nx - 1);
}

/**
 * @brief Perform one Crank--Nicolson time step. 
 * 
 * 1. Assemble RHS = (I + r/2 A) u^n.
 * 2. Add contributions from Dirichlet boundary.
 * 3. Solve tridiagonal system (I - r/2 A) u^{n+1} = RHS.
 * 4. Apply boundary conditions to u.
 * 5. Copy updated interior solution back into u. 
 * 
 * @param t Current time
 */
void CrankNicolson1D::step(double t)
{
    // Step 1: Build RHS using u^n
    for (std::size_t i = 1; i < Nx; i++) 
    {
        rhs[i - 1] = (1.0 - r) * u[i] + (r / 2.0) * (u[i - 1] + u[i + 1]);
    }

    // Step 2: Add Dirichlet boundary contributions
    rhs[0]    += (r / 2.0) * bc.left(t);
    rhs[Nx-2] += (r / 2.0) * bc.right(t);

    // Step 3: Solve tridiagonal system
    tridiagonalSolve(lower, diag, upper, rhs);

    // Step 4: Enforce Dirichlet boundary conditions
    applyBoundaryConditions(t);

    // Step 5: Copy interior solution back into u. 
    for (std::size_t i = 1; i < Nx; i++) 
    {
        u[i] = rhs[i - 1];
    }
}