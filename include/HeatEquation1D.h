#pragma once

/**
 * @brief Represents the 1D heat equation ∂u/∂t = α ∂²u/∂x²
 * 
 * This class stores and provides access to the diffusion constant α. 
 */
class HeatEquation1D {
public: 
    /**
     * @brief Construct a HeatEquation1D object.
     *
     * @param alpha Diffusion coefficient (α > 0)
     */
    HeatEquation1D(double alpha);
    
    /**
     * @brief Return the diffusion coefficient α.
     * @return α
     */
    double diffusion() const;

private:
    double alpha; //< Diffusion coefficient
};