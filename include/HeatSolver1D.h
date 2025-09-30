#pragma once

#include "Grid1D.h"
#include "HeatEquation1D.h"
#include "BoundaryCondition.h"
#include "InitialCondition.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

/**
 * @brief Abstract base class for 1D heat equation solvers.
 * 
 * This class provides a common interface for time-stepping solvers of the 1D heat equation.
 * Specific methods (Forward Euler, Backward Euler, Crank-Nicolson) should inherit
 * from this class and implement the pure virtual step() method.
 */
class HeatSolver1D {
public:
    /**
     * @brief Construct a HeatSolver1D object. 
     * 
     * @param grid Reference to the spatial grid.
     * @param pde Reference to the HeatEquation1D object (stores diffusion a)
     * @param bc Reference to the boundary condition.
     * @param ic Reference to the initial condition.
     * @param dt Time step size.
     * @param T Final simulation time. 
     * 
     * @note This constructor does not initialize the solution vector. 
     *       Call initialize() before running the simulation. 
     */
    HeatSolver1D(const Grid1D& grid,
                 const HeatEquation1D& pde,
                 const BoundaryCondition& bc,
                 const InitialCondition& ic,
                 double dt, 
                 double T);
    
    /**
     * @brief Initialize the solver.
     * 
     * Sets up the solution vector. 
     */
    virtual void initialize();

    /**
     * @brief Run the simulation from t=0 to t=T.
     */
    virtual void runSimulation();
    
    /**
     * @brief Run the simulation and save solution snapshots periodically.
     * 
     * Saves the solution vector at specified intervals during the time integration.
     * 
     * @param prefix Prefix for output filenames (e.g., "heat_FE").
     * @param saveEvery Save every N timesteps (default = 10). 
     */
    void runSimulationWithSnapshots(const std::string& prefix, int saveEvery = 10);

    /**
     * @brief Get the current solution vector.
     * 
     * @return const reference to vector of solution values u. 
     */
    const std::vector<double>& getSolution() const;

    /**
     * @brief Save the current solution to a file.
     * 
     * @param filename Name of output file. 
     */
    void save(const std::string& filename) const; 

protected:
    const Grid1D& grid;                 //< Spatial grid
    const HeatEquation1D& pde;          //< PDE parameters (α)
    const BoundaryCondition& bc;        //< Boundary conditions
    const InitialCondition& ic;         //< Initial conditions

    std::size_t Nx, Nt;                 //< Number of spatial and time points.
    double dt, T, dx;                   //< Time step, final time, grid spacing

    std::vector<double> u;              //< Solution vector at current time step. 

    /**
     * @brief Perform a single time step (pure virtual).
     * 
     * Must be implemented by derived classes to update u.
     * 
     * @param t Current time.  
     */
    virtual void step(double t) = 0;

    /**
     * @brief Apply boundary conditions to a given vector. 
     * 
     * @param t Current time.  
     */
    void applyBoundaryConditions(double t);

    /**
     * @brief Apply boundary conditions to a given vector.
     *
     * @param t Current time
     * @param u_vec Vector to which boundary conditions are applied
     */
    void applyBoundaryConditionsUpdated(double t, std::vector<double>& u_vec);
};