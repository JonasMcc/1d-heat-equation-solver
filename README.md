# A modular C++ framework for solving the 1D heat equation
This project provides a flexible C++ framework for numerically solving the 1D heat equation
\[
    \frac{\partal u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2}
\]
where $x \in [a, b]$, $t \geq 0$, and $\alpha$ is the diffusion coefficient. The following parameters are fully configurable:
- Spatial domain and grid resolution. 
- Dirichlet boundary conditions. 
- Initial condition.
- Diffusion coefficient $\alpha$.

The framework includes three finite difference methods: 
- Forward Euler
- Backward Euler
- Crank–Nicolson

The Backward Euler and Crank–Nicolson methods are solved using a tridiagonal solver implementing the Thomas algorithm. 

## Example




