#include "InitialCondition.h"

/**
 * @brief Construct a InitialCondition object. 
 */
InitialCondition::InitialCondition(std::function<double(double)> f_) : f(f_) {}

/**
 * @brief Evaluate the initial function at a spatial coordinate x. 
 */
double InitialCondition::operator()(double x) const { return f(x); }