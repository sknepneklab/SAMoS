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
  p.Nx = Nx; p.Ny = Ny; p.Nz = Nz;
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

