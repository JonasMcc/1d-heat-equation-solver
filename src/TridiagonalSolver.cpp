#include "TridiagonalSolver.h"

/**
 * @brief Solve a tridiagonal linear system Ax = d using the Thomas algorithm.
 * 
 * The algorithm consists of: 
 * 1. Forward elimination: reduce the system to upper-triangular form.
 * 2. Back substitution: solve for the unknown once upper-triangular form has been achieved.
 * 
 * Complexity: O(n)
 * 
 * @param a Sub-diagonal coefficients (length n-1)
 * @param b Main diagonal coefficients (length n)
 * @param c Super-diagonal coefficients (length n-1)
 * @param d Right-hand side vector (length n) 
 */
void tridiagonalSolve(const std::vector<double>& a,
                      const std::vector<double>& b,
                      const std::vector<double>& c, 
                      std::vector<double>& d)
{
    std::size_t n = d.size();
    std::vector<double> c_star(n-1), d_star(n);

    // Step 1: Forward elimination
    // Compute modified coefficients c* and d* such that system becomes upper-triangular
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (std::size_t i = 1; i < n; i++) {
        // Elimination factor
        double m = 1.0 / (b[i] - a[i-1]*c_star[i-1]);

        // Update c* (only until second-to-last equation)
        if (i < n-1) 
        {
            c_star[i] = c[i]*m; 
        }
        
        // Update d*
        d_star[i] = (d[i] - a[i-1]*d_star[i-1])*m;
    }

    // Step 2: Back substitution
    // Solve the last variable directly
    d[n-1] = d_star[n-1];

    // Recurse backwards to solve for all variables
    for (int i = static_cast<int>(n) - 2; i >=0; i--)
    {
         d[i] = d_star[i] - c_star[i]*d[i+1];
    }
}