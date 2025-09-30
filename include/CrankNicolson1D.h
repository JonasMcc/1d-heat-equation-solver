#pragma once
#include "HeatSolver1D.h"
#include "TridiagonalSolver.h"

/**
 * @brief Crank--Nicolson solver for the 1D heat equation.
 * 
 * Implements the semi-implicit Crank--Nicolson method:
 * 
 *      (I - r/2 A) u^{n+1} = (I + r/2 A) u^n,
 * 
 * where r = α Δt / Δx² and A is the standard second-derivative finite-difference matrix.
 * 
 * @note Unlike Forward Euler, Crank--Nicolson is unconditionally stable and second-order accurate in time. 
 */
class CrankNicolson1D : public HeatSolver1D {
public:
    /**
     * @brief Construct a CrankNicolson1D solver.
     * 
     * Precomputes the left-hand side tridiagonal coefficients (I - r/2 A).
     * 
     * @param grid Spatial grid
     * @param pde Heat equation parameters (diffusion)
     * @param bc Boundary condition
     * @param ic Initial condition
     * @param dt Time step size
     * @param T Final simulation time
     */
    CrankNicolson1D(const Grid1D& grid,
                    const HeatEquation1D& pde,
                    const BoundaryCondition& bc,
                    const InitialCondition& ic,
                    double dt,
                    double T);
protected:
    /**
     * @brief Perform a single Crank--Nicolson time step.
     * 
     * This function does the following:
     *  - Builds the right-hand side (I + r/2 A) u^n.
     *  - Applies the boundary conditions
     *  - Solves the tridiagonal system.
     *  - Updates the solution
     * 
     * @param t Current time
     */
    void step(double t) override; 

private:
    double r;                                   //< CFL number
    std::vector<double> lower, diag, upper;     //< Tridiagonal of LHS
    std::vector<double> rhs;                    //< RHS vector
};