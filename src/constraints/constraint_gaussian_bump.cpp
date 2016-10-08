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
  p.Nx = Nx; p.Ny = Ny; p.Nz = Nz;
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

