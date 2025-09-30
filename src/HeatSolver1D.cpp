#include "HeatSolver1D.h"
#include <filesystem>

/**
 * @brief Construct a HeatSolver1D object.
 *
 * Initializes grid, PDE, boundary conditions, initial condition,
 * time step, and final time. Allocates the solution vector u.
 *
 * @param grid_ Spatial grid
 * @param pde_ Heat equation parameters (diffusion coefficient α)
 * @param bc_ Boundary conditions
 * @param ic_ Initial condition
 * @param dt_ Time step size
 * @param T_  Final simulation time
 */
HeatSolver1D::HeatSolver1D(
    const Grid1D& grid_, 
    const HeatEquation1D& pde_, 
    const BoundaryCondition& bc_,
    const InitialCondition& ic_,
    double dt_,
    double T_)
    : grid(grid_)
    , pde(pde_)
    , bc(bc_)
    , ic(ic_)
    , dt(dt_)
    , T(T_)
{
    Nx = grid.size() - 1;                       //< Number of intervals
    dx = grid.spacing();                        //< Spacing between grid points
    Nt = static_cast<std::size_t>(T / dt);      //< Number of time steps

    u.resize(grid.size(), 0.0);                 //< Initialize solution vector.
}

/**
 * @brief Apply initial confition to the solution vector. 
 * 
 * This sets u(x, 0) = ic(x) for all spatial grid points.
 */
void HeatSolver1D::initialize() 
{
    const auto& x = grid.points();

    for (std::size_t i = 0; i <= Nx; i++)
    {
        u[i] = ic(x[i]);
    }
}

/**
 * @brief Run the time-step simulation. 
 * 
 * Calls the derived class step() function for each time step.
 * Boundary conditions are expected to be applied inside step() or separately.
 */
void HeatSolver1D::runSimulation()
{
    for(std::size_t n = 0; n < Nt; n++)
    {
        double t = n * dt;
        step(t);    // Pure virtual function implemented by derived solver
    }
}

/**
 * @brief Run the simulation and periodically save snapshots of the solution.
 *
 * Generates output files named like "<prefix>_stepN.dat" at the chosen interval.
 *
 * @param prefix Prefix for filenames
 * @param saveEvery Save every N time steps
 */
void HeatSolver1D::runSimulationWithSnapshots(const std::string& prefix, int saveEvery)
{
    for (std::size_t n = 0; n <= Nt; ++n)
    {
        double t = n * dt;

        // Save snapshot
        if (n % saveEvery == 0)
        {
            std::string filename = prefix + "_step" + std::to_string(n) + ".dat";
            save(filename);
        }

        // Skip stepping after final snapshot
        if (n < Nt)
        {
            step(t);
        }
    }
}

/**
 * @brief Save the solution vector to a file.
 * 
 * Output format: x_i u_i
 * 
 * @param filename Filename of output file. 
 */
void HeatSolver1D::save(const std::string& filename) const 
{
    std::filesystem::create_directories("results");

    std::ofstream out("results/" + filename);

    for (std::size_t i = 0; i < grid.size(); i++)
    {
        out << grid[i] << " " << u[i] << "\n";
    }
}

//void HeatSolver1D::save(const std::string& filename) const 
//{
//    std::ofstream fout(filename);
//    const auto& x = grid.points();

//    for (std::size_t i = 0; i <= Nx; i++) 
    //{
      //  fout << x[i] << " " << u[i] << "\n";
    //}

    //fout.close();
    //std::cout << "Results saved to " << filename << "\n";
//}

/**
 * @brief Apply Dirichlet boundary conditions to the solution vector u.
 *
 * Sets u[0] and u[Nx] to the left and right boundary values.
 *
 * @param t Current time
 */
void HeatSolver1D::applyBoundaryConditions(double t)
{
    u[0]  = bc.left(t);
    u[Nx] = bc.right(t);
}

/**
 * @brief Apply boundary conditions to a provided vector.
 *
 * Useful if intermediate vectors are used in solver implementations.
 *
 * @param t Current time
 * @param u_vec Vector to which boundary conditions are applied
 */
void HeatSolver1D::applyBoundaryConditionsUpdated(double t, std::vector<double>& u_vec)
{
    u_vec[0]  = bc.left(t);
    u_vec[Nx] = bc.right(t);
}