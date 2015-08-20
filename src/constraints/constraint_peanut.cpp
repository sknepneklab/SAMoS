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
 * \file constraint_peanut.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 11-Nov-2014
 * \brief Implementation of the peanut constraint
 */ 

#include "constraint_peanut.hpp"


/*! Compute normal to the surface at p
*/
void ConstraintPeanut::compute_normal(Particle& p, double& Nx, double& Ny, double& Nz)
{
  this->compute_gradient(p,Nx,Ny,Nz);
  // Normalize N
  double len_N = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
  Nx /= len_N; Ny /= len_N; Nz /= len_N;
}

/*! Compute gradient at a point 
 *  \param p reference to a point
*/
void ConstraintPeanut::compute_gradient(Particle&p, double& gx, double& gy, double& gz)
{
  double x = p.x, y = p.y, z = p.z;
  double fact_1 = 4.0*m_a*m_a;
  double fact_2 = 4.0*(x*x + y*y + z*z);
  gx = (-fact_1 + fact_2)*x;
  gy = ( fact_1 + fact_2)*y; 
  gz = ( fact_1 + fact_2)*z;
}

/*! Compute constraint value at particle p
 *  \param p reference to the particle
*/
double ConstraintPeanut::constraint_value(Particle& p)
{
  double x = p.x, y = p.y, z = p.z;
  double tr = x*x + y*y + z*z;
  double g = m_a*m_a*m_a*m_a - m_b*m_b*m_b*m_b - 2.0*m_a*m_a*(x*x - y*y - z*z) + tr*tr;
  return g;
}


