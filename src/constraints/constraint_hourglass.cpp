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
 *   (c) 2013, 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file constraint_hourglass.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 09-Feb-2015
 * \brief Implementation of the hourglass constraint
 */ 

#include "constraint_hourglass.hpp"

/*! Compute normal to the surface at p
*/
void ConstraintHourglass::compute_normal(Particle& p, double& Nx, double& Ny, double& Nz)
{
  // Compute unit gradient of the implicit function (normal will be normalized gradient)
  this->compute_gradient(p,Nx,Ny,Nz);
  // Normalize N
  double len_N = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
  Nx /= len_N; Ny /= len_N; Nz /= len_N;
}

/*! Compute gradient at a point 
 *  \param p reference to a point
*/
void ConstraintHourglass::compute_gradient(Particle&p, double& gx, double& gy, double& gz)
{
  double x = p.x, y = p.y, z = p.z;
  BoxPtr box = m_system->get_box();
  double lz = box->Lz;
  gx = 2.0*x;
  gy = 2.0*y; 
  gz = (-4*m_A*m_n*M_PI*cos((2*m_n*M_PI*z)/lz)*(m_R + m_A*sin((2*m_n*M_PI*z)/lz)))/lz;
}

/*! Compute constraint value at particle p
 *  \param p reference to the particle
*/
double ConstraintHourglass::constraint_value(Particle& p)
{
  double x = p.x, y = p.y, z = p.z;
  BoxPtr box = m_system->get_box();
  double lz = box->Lz;
  double zfact = m_R + m_A*sin(2.0*M_PI/lz*m_n*z);
  double g = x*x + y*y - zfact*zfact;
  return g;
}

