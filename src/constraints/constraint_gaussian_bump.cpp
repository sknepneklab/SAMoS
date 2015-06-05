/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014, 2015
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file constraint_gaussian_bump.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Jun-2015
 * \brief Implementation of the toroidal constraint
 */ 

#include "constraint_gaussian_bump.hpp"


/*! Compute normal to the surface at p
 *  
 *  Normal vector is computed as the normalized gradient vector at point p.
 *  Coordinates of gradient vector are given as:
 * 
 *  \f$ g_x = -\frac{2 A x e^{-\frac{x^2}{a^2}-\frac{y^2}{b^2}}}{a^2}, \f$
 *  \f$ g_y = -\frac{2 A y e^{-\frac{x^2}{a^2}-\frac{y^2}{b^2}}}{b^2}, \f$
 *  \f$ g_z = -1. \f$  
 *  
 *  \param p Reference to the particle object
 *  \param Nx x coordinate of the normal
 *  \param Ny y coordinate of the normal
 *  \param Nz z coordinate of the normal
*/
void ConstraintGaussianBump::compute_normal(Particle& p, double& Nx, double& Ny, double& Nz)
{
  this->compute_gradient(p,Nx,Ny,Nz);
  double len_N = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
  Nx /= len_N; Ny /= len_N; Nz /= len_N;
}

/*! Compute gradient at a point 
 *  \param p reference to a point
*/
void ConstraintGaussianBump::compute_gradient(Particle&p, double& gx, double& gy, double& gz)
{
  double x = p.x, y = p.y;
  double fact = exp(-x*x/m_a-y*y/m_b); 
  gx = -2.0*m_A*x*fact/m_a;
  gy = -2.0*m_A*y*fact/m_b;
  gz = -1.0;
}

/*! Compute constraint value at particle p
 *  \param p reference to the particle
*/
double ConstraintGaussianBump::constraint_value(Particle& p)
{
  double x = p.x, y = p.y, z = p.z;
  double g = m_A*exp(-x*x/m_a - y*y/m_b)-z;  // implicit equation for gaussian_bump
  return g;
}

