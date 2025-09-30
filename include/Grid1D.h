#pragma once
#include <vector>

/**
 * @brief Represents a uniform 1D grid on the interval [a, b].
 * ADD
 */
class Grid1D {
public:
    /**
     * @brief Construct a 1D grid on the interval [a, b].
     *
     * @param N Number of intervals. The grid will have N+1 points.
     * @param a Left boundary of the domain.
     * @param b Right boundary of the domain.
     *
     * @note The grid spacing is computed as dx = (b - a) / N.
     *       The points are distributed uniformly from a to b.
     */
    Grid1D(int N, double a, double b);

    /**
     * @brief Get the number of grid points.
     * @return Total number of points in the grid (N+1).
     */
    std::size_t size() const;

    /**
     * @brief Get the grid spacing.
     * @return dx = (b - a) / N
     */
    double spacing() const;

    /**
     * @brief Get the left boundary of the domain.
     * @return Value of a
     */
    double xmin() const;

    /**
     * @brief Get the right boundary of the domain.
     * @return Value of b
     */
    double xmax() const;

    /**
     * @brief Get a const reference to the grid points.
     * @return Vector of grid coordinates.
     */
    const std::vector<double>& points() const; 

    /**
     * @brief Access grid point by index.
     *
     * @param i Index of the grid point (0 ≤ i < size()).
     * @return Coordinate value x[i].
     */
    double operator[](std::size_t i) const; 

private:
    int N;                 //< Number of intervals
    double a, b, dx;       //< Domain boundaries and spacing
    std::vector<double> x; //< Vector of coordinates of grid points
};