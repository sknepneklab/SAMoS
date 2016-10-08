/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

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
  p.Nx = Nx; p.Ny = Ny; p.Nz = Nz;
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


