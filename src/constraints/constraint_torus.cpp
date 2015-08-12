/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file constraint_torus.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 15-Dec-2014
 * \brief Implementation of the toroidal constraint
 */ 

#include "constraint_torus.hpp"


/*! Compute normal to the surface at p
 *  
 *  Normal vector is computed as the normalized gradient vector at point p.
 *  Coordinates of gradient vector are given as:
 * 
 *  \f$ g_x = -\frac{2x\left(c-\sqrt{x^2+y^2}\right)}{\sqrt{x^2+y^2}}, \f$
 *  \f$ g_y = -\frac{2y\left(c-\sqrt{x^2+y^2}\right)}{\sqrt{x^2+y^2}}, \f$
 *  \f$ g_z = 2z. \f$  
 *  
 *  \param p Reference to the particle object
 *  \param Nx x coordinate of the normal
 *  \param Ny y coordinate of the normal
 *  \param Nz z coordinate of the normal
*/
void ConstraintTorus::compute_normal(Particle& p, double& Nx, double& Ny, double& Nz)
{
  this->compute_gradient(p,Nx,Ny,Nz);
  double len_N = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
  Nx /= len_N; Ny /= len_N; Nz /= len_N;
}

/*! Compute gradient at a point 
 *  \param p reference to a point
*/
void ConstraintTorus::compute_gradient(Particle&p, double& gx, double& gy, double& gz)
{
  double x = p.x, y = p.y, z = p.z;
  double sqr_xy = sqrt(x*x + y*y);
  double fact = -2.0*(m_c - sqr_xy)/sqr_xy; 
  gx = x*fact;
  gy = y*fact;
  gz = 2.0*z;
}

/*! Compute constraint value at particle p
 *  \param p reference to the particle
*/
double ConstraintTorus::constraint_value(Particle& p)
{
  double x = p.x, y = p.y, z = p.z;
  double ft = m_c - sqrt(x*x + y*y);
  double g = ft*ft + z*z - m_a*m_a;  // implicit equation for torus
  return g;
}


