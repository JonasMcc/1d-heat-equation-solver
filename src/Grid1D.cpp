#include "Grid1D.h"

/**
 * @brief Constructs a uniform 1D grid. 
 * 
 * Fills the vector of grid points from a to b with N+1 points.
 * 
 * @throws std::invalid_argument if N <= 0 and if b is not > a. 
 */
Grid1D::Grid1D(
    int N_,         
    double a_, 
    double b_)
    : N(N_)
    , a(a_)
    , b(b_)
    , dx((b_ - a_)/ N_)
    , x(N_ + 1)
{
    if (N_ <= 0)
    {
        throw std::invalid_argument("Grid1D: N must be > 0.");
    }

    if (b_ <= a_)
    {
        throw std::invalid_argument("Grid1D: b must be > a.");
    }

    for (std::size_t i = 0; i <= static_cast<std::size_t>(N); i++)
    {
        x[i] = a + i * dx;
    }
}

/**
 * @brief Get the number of grid points.
 * @return The total number of grid points (N+1). 
 */
std::size_t Grid1D::size() const { return static_cast<std::size_t>(N) + 1; }

/**
* @brief Get the spacing dx between the grid points. 
* @return dx = (b - a) / N
*/
double Grid1D::spacing() const { return dx; }
    
/**
* @brief Get the left boundary of the domain.
* @return Value of a
*/
double Grid1D::xmin() const { return a; }

/**
* @brief Get the right boundary of the domain.
* @return Value of b
*/
double Grid1D::xmax() const { return b; }

/**
* @brief Get a const reference to the grid points.
* @return Vector of grid coordinates.
*/
const std::vector<double>& Grid1D::points() const { return x; }

/**
 * @brief Access grid point by index. 
 */
double Grid1D::operator[](std::size_t i) const { return x[i]; }