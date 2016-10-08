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
  p.Nx = Nx; p.Ny = Ny; p.Nz = Nz;
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

