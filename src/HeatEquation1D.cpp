#include "HeatEquation1D.h"
#include <stdexcept>

/**
 * @brief Construct a HeatEquation1D object.
 * 
 * This class stores the diffusion coefficient for later use by solvers.
 * 
 * @throws std::invalid_argument if diffusion constant is not > 0. 
 */
HeatEquation1D::HeatEquation1D(double alpha_) : alpha(alpha_)
{
    if (!(alpha_ > 0.0)) 
    {
        throw std::invalid_argument("HeatEquation1D: The diffusion constant must be > 0.");
    }
}

/**
 * @brief Return the diffusion coefficient. 
 */
double HeatEquation1D::diffusion() const { return alpha; }
