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
 * \file constraint_tetrahedron.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Aug-2015
 * \brief Implementation of the tetrahedron constraint
 */ 

#include "constraint_tetrahedron.hpp"

/*! Compute normal to the surface at p
*/
void ConstraintTetrahedron::compute_normal(Particle& p, double& Nx, double& Ny, double& Nz)
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
void ConstraintTetrahedron::compute_gradient(Particle&p, double& gx, double& gy, double& gz)
{
  Vector3d pos = Vector3d(p.x,p.y,p.z);
  Vector3d p_a1 = pos - m_a1;
  Vector3d p_a2 = pos - m_a2;
  Vector3d p_a3 = pos - m_a3;
  Vector3d p_a4 = pos - m_a4;
  double l1 = p_a1.len();
  double l2 = p_a2.len();
  double l3 = p_a3.len();
  double l4 = p_a4.len();
  double l1_3 = l1*l1*l1, l2_3 = l2*l2*l2, l3_3 = l3*l3*l3, l4_3 = l4*l4*l4;
  gx = -(p_a1.x/l1_3 + p_a2.x/l2_3 + p_a3.x/l3_3 + p_a4.x/l4_3);
  gy = -(p_a1.y/l1_3 + p_a2.y/l2_3 + p_a3.y/l3_3 + p_a4.y/l4_3);
  gz = -(p_a1.z/l1_3 + p_a2.z/l2_3 + p_a3.z/l3_3 + p_a4.z/l4_3);
}

/*! Compute constraint value at particle p
 *  The constraint is given as the equipotential surface of the potential
 *  \f$ V(\vec r) = \sum_{i=1}^{4}\frac{1}{|\vec r - \vec a_i|} \f$, where 
 *  \f$ \vec a_i \f$ are position of the four unit charges generating the 
 *  potential.
 *  \param p reference to the particle
*/
double ConstraintTetrahedron::constraint_value(Particle& p)
{
  Vector3d pos = Vector3d(p.x,p.y,p.z);
  Vector3d p_a1 = pos - m_a1;
  Vector3d p_a2 = pos - m_a2;
  Vector3d p_a3 = pos - m_a3;
  Vector3d p_a4 = pos - m_a4;
  double g = 1.0/sqrt(dot(p_a1, p_a1));
  g += 1.0/sqrt(dot(p_a2, p_a2));
  g += 1.0/sqrt(dot(p_a3, p_a3));
  g += 1.0/sqrt(dot(p_a4, p_a4));
  return g - m_value;
}

/*! Rescale tetrahedron size and make sure that all particles are still on it.
 *  Rescaling is done only at certain steps and only if rescale 
 *  factor is not equal to 1.
 */
bool ConstraintTetrahedron::rescale()
{
  if (m_rescale != 1.0)
  {
    int step = m_system->get_step();
    if ((step % m_rescale_freq == 0) && (step < m_rescale_steps))
    {
      m_s *= m_scale;
      m_a1 = m_s*tetra_a1;
      m_a2 = m_s*tetra_a2;
      m_a3 = m_s*tetra_a3;
      m_a4 = m_s*tetra_a4;
      for  (int i = 0; i < m_system->size(); i++)
      {
        Particle& p = m_system->get_particle(i);
        this->enforce(p);
      }
      return true;
    }
  }
  return false;
}

