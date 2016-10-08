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
  p.Nx = Nx; p.Ny = Ny; p.Nz = Nz;
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


