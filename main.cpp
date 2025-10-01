#include "Grid1D.h"
#include "HeatEquation1D.h"
#include "BoundaryCondition.h"
#include "InitialCondition.h"
#include "ForwardEuler1D.h"
#include "BackwardEuler1D.h"
#include "CrankNicolson1D.h"
#include <iostream>
#include <cmath>

int main()
{
    // Discretization 
    int Nx = 50;
    double a = 0.0;
    double b = 1.0;
    double alpha = 1.0;
    double dt = 0.0001;
    double T = 0.1;

    // Grid
    Grid1D grid(Nx, a, b);

    // Heat equation with diffusion alpha. 
    HeatEquation1D pde(alpha);

    // Boundary conditions: u(0,t)=0, u(1,t)=0
    BoundaryCondition bc(
        [](double /*t*/) { return 0.0; },
        [](double /*t*/) { return 0.0; }
    );

    // Boundary condition with sine wave.
    BoundaryCondition bc_sine(
        [](double t) { return std::sin(t); },
        [](double /*t*/) { return 0.0; }
    );

    // Initial condition: sine wave
    InitialCondition ic([](double x) {
        return std::sin(M_PI * x);
    });

    // Choose a scheme. We will use Crank--Nicolson. 
    ForwardEuler1D solverFE(grid, pde, bc, ic, dt, T);
    BackwardEuler1D solverBE(grid, pde, bc, ic, dt, T);
    CrankNicolson1D solverCN(grid, pde, bc, ic, dt, T);

    // With snapshots
    solverFE.initialize();
    solverFE.runSimulationWithSnapshots("HeatFE", 50);

    // Initialize and run FE
    //solverFE.initialize();
    //solverFE.runSimulation();
    //solverFE.save("HeatFE.dat");

    // Initialize and run CN
    //solverBE.initialize();
    //solverBE.runSimulation();
    //solverBE.save("HeatBE.dat");

    // Initialize and run CN
    //solverCN.initialize();
    //solverCN.runSimulation();
    //solverCN.save("HeatCN.dat");

    std::cout << "Simulation finished. Results in ./results/ directory.\n";
    
    return 0;
}