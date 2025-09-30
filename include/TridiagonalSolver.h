#pragma once
#include <vector>

/**
 * @brief Solve a tridiagonal linear system Ax = d using the Thomas algorithm. 
 * 
 * The tridiagonal system has the form:
 *
 *   b[i] * x[i] + c[i] * x[i+1] + a[i-1] * x[i-1] = d[i],   for i = 0,...,n-1
 *
 * where:
 * - a (sub-diagonal) has length n-1
 * - b (main diagonal) has length n
 * - c (super-diagonal) has length n-1
 * - d (right-hand side) has length n
 * 
 * @param a Sub-diagonal coefficients (length n-1)
 * @param b Main diagonal coefficients (length n)
 * @param c Super-diagonal coefficients (length n-1)
 * @param rhs Right-hand side vector (length n) 
 * 
 * @note This function overrides rhs with the solution.  
 */
void tridiagonalSolve(const std::vector<double>& a,
                      const std::vector<double>& b,
                      const std::vector<double>& c, 
                      std::vector<double>& rhs);