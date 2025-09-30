#include "HeatEquation1D.h"

/**
 * @brief Construct a HeatEquation1D object.
 * 
 * This class stores the diffusion coefficient for later use by solvers.
 */
HeatEquation1D::HeatEquation1D(double alpha_) : alpha(alpha_) {}

/**
 * @brief Return the diffusion coefficient. 
 */
double HeatEquation1D::diffusion() const { return alpha; }
