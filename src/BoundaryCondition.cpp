#include "BoundaryCondition.h"

/**
 * @brief Construct a BoundaryCondition object.
 *
 * Stores the functions for the left and right boundaries.
 */
BoundaryCondition::BoundaryCondition(std::function<double(double)> leftFunction_,
                                     std::function<double(double)> rightFunction_)
    : leftFunction(leftFunction_)
    , rightFunction(rightFunction_)
{
}

/**
 * @brief Evaluate the left boundary function at time t.
 */
double BoundaryCondition::left(double t) const { return leftFunction(t); }

/**
 * @brief Evaluate the right boundary function at time t.
 */
double BoundaryCondition::right(double t) const { return rightFunction(t); }