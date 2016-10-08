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
 * \file constraint_ellipsoid.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 15-Dec-2014
 * \brief Implementation of the toroidal constraint
 */ 

#include "constraint_ellipsoid.hpp"


/*! Compute normal to the surface at p
 *  
 *  Normal vector is computed as the normalized gradient vector at point p.
 *  Coordinates of gradient vector are given as:
 * 
 *  \f$ g_x = \frac{2x}{a^2}, \f$
 *  \f$ g_y = \frac{2y}{b^2}, \f$
 *  \f$ g_z = \frac{2z}{c^2}. \f$  
 *  
 *  \param p Reference to the particle object
 *  \param Nx x coordinate of the normal
 *  \param Ny y coordinate of the normal
 *  \param Nz z coordinate of the normal
*/
void ConstraintEllipsoid::compute_normal(Particle& p, double& Nx, double& Ny, double& Nz)
{
  this->compute_gradient(p,Nx,Ny,Nz);
  double len_N = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
  Nx /= len_N;  Ny /= len_N;  Nz /= len_N;
  p.Nx = Nx; p.Ny = Ny; p.Nz = Nz;
}

/*! Computes gradient onto the surface
 *  \param p Reference to the particle object
 *  \param gx x component of the gradient vector (returned)
 *  \param gy y component of the gradient vector (returned)
 *  \param gz z component of the gradient vector (returned)
*/
void ConstraintEllipsoid::compute_gradient(Particle& p, double& gx, double& gy, double& gz)
{
  gx = 2.0*p.x/m_a/m_a;
  gy = 2.0*p.y/m_b/m_b;
  gz = 2.0*p.z/m_c/m_c;
}


/*! Computes value of the ellipsoid constraint for particle 
 *  \param p reference to a particle
*/
double ConstraintEllipsoid::constraint_value(Particle& p)
{
  double x = p.x, y = p.y, z = p.z;
  return x*x/m_a/m_a + y*y/m_b/m_b + z*z/m_c/m_c - 1.0;
}

