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
 * \file constraint_gyroid.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 15-Dec-2014
 * \brief Implementation of the toroidal constraint
 */ 

#include "constraint_gyroid.hpp"

/*! Compute normal to the surface at p
 *  
 *  Normal vector is computed as the normalized gradient vector at point p.
 *  Coordinates of gradient vector are given as:
 * 
 *  \f$ g_x = \frac{2 \pi  \left(\cos \left(\frac{2 \pi  x}{lx}\right) \cos \left(\frac{2 \pi  z}{lz}\right)-\sin \left(\frac{2 \pi  x}{lx}\right) \sin \left(\frac{2 \pi 
   y}{ly}\right)\right)}{lx}, \f$
 *  \f$ g_y = \frac{2 \pi  \left(\cos \left(\frac{2 \pi  x}{lx}\right) \cos \left(\frac{2 \pi  y}{ly}\right)-\sin \left(\frac{2 \pi  y}{ly}\right) \sin \left(\frac{2 \pi 
   z}{lz}\right)\right)}{ly}, \f$
 *  \f$ g_z = \frac{2 \pi  \left(\cos \left(\frac{2 \pi  y}{ly}\right) \cos \left(\frac{2 \pi  z}{lz}\right)-\sin \left(\frac{2 \pi  x}{lx}\right) \sin \left(\frac{2 \pi 
   z}{lz}\right)\right)}{lz}. \f$  
 *  
 *  \param p Reference to the particle object
 *  \param Nx x coordinate of the normal
 *  \param Ny y coordinate of the normal
 *  \param Nz z coordinate of the normal
*/
void ConstraintGyroid::compute_normal(Particle& p, double& Nx, double& Ny, double& Nz)
{
  this->compute_gradient(p,Nx,Ny,Nz);
  double len_N = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
  Nx /= len_N; Ny /= len_N; Nz /= len_N;
}


/*! Compute gradient at a point 
 *  \param p reference to a point
*/
void ConstraintGyroid::compute_gradient(Particle&p, double& gx, double& gy, double& gz)
{
  double x = p.x, y = p.y, z = p.z;
  BoxPtr box = m_system->get_box();
  double lx = box->Lx, ly = box->Ly, lz = box->Lz;
  gx = (2*M_PI*(cos((2*M_PI*x)/lx)*cos((2*M_PI*z)/lz) - sin((2*M_PI*x)/lx)*sin((2*M_PI*y)/ly)))/lx;
  gy = (2*M_PI*(cos((2*M_PI*x)/lx)*cos((2*M_PI*y)/ly) - sin((2*M_PI*y)/ly)*sin((2*M_PI*z)/lz)))/ly;
  gz = (2*M_PI*(cos((2*M_PI*y)/ly)*cos((2*M_PI*z)/lz) - sin((2*M_PI*x)/lx)*sin((2*M_PI*z)/lz)))/lz;
}

/*! Compute constraint value at particle p
 *  \param p reference to the particle
*/
double ConstraintGyroid::constraint_value(Particle& p)
{
  double x = p.x, y = p.y, z = p.z;
  BoxPtr box = m_system->get_box();
  double lx = box->Lx, ly = box->Ly, lz = box->Lz;
  return cos((2*M_PI*z)/lz)*sin((2*M_PI*x)/lx) + cos((2*M_PI*x)/lx)*sin((2*M_PI*y)/ly) + cos((2*M_PI*y)/ly)*sin((2*M_PI*z)/lz);
}

