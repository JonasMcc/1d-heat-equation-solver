#include "TridiagonalSolver.h"
#include <stdexcept>
#include <limits>
#include <cmath>

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
 * @param rhs Right-hand side vector (length n) 
 * 
 * @throws std::invalid_argument if input vector sizes are inconsistent.
 * @throws std::runtime_error.   if a zero pivot is encountered (singular system).
 */
void tridiagonalSolve(const std::vector<double>& a,
                      const std::vector<double>& b,
                      const std::vector<double>& c, 
                      std::vector<double>& rhs)
{
    std::size_t n = rhs.size();

    // Basic size checks
    if (n == 0) 
    {
        throw std::invalid_argument("tridiagonalSolve: system size n == 0");
    }
    if (b.size() != n) 
    {
        throw std::invalid_argument("tridiagonalSolve: b.size() != d.size()");
    }
    if (n > 1) 
    {
        if (a.size() != n - 1 || c.size() != n - 1) 
        {
            throw std::invalid_argument("tridiagonalSolve: a and c must have length n-1 when n>1");
        }
    } 
    else 
    {
        // n == 1, a and c should be empty
        if (!a.empty() || !c.empty()) 
        {
            throw std::invalid_argument("tridiagonalSolve: a/c must be empty for n==1");
        }
    }

    // Temporary vectors
    std::vector<double> c_star(n - 1);
    std::vector<double> d_star(n);

    // First equation
    if (std::abs(b[0]) < std::numeric_limits<double>::epsilon()) 
    {
        throw std::runtime_error("tridiagonalSolve: zero pivot in b[0]");
    }
    
    c_star[0] = (n > 1) ? c[0] / b[0] : 0.0;
    d_star[0] = rhs[0] / b[0];

    // Forward elimination
    for (std::size_t i = 1; i < n; ++i) 
    {
        double denom = b[i] - a[i - 1] * c_star[i - 1];

        if (std::abs(denom) < std::numeric_limits<double>::epsilon()) 
        {
            throw std::runtime_error("tridiagonalSolve: zero pivot encountered during elimination");
        }
        
        double m = 1.0 / denom;

        if (i < n - 1) 
        {
            c_star[i] = c[i] * m;
        }
        
        d_star[i] = (rhs[i] - a[i - 1] * d_star[i - 1]) * m;
    }

    // Back substitution: overwrite d with the solution
    rhs[n - 1] = d_star[n - 1];

    for (int i = static_cast<int>(n) - 2; i >= 0; --i) 
    {
        rhs[i] = d_star[i] - c_star[i] * rhs[i + 1];
    }
}