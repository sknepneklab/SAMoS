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
 * \file constraint_cylinder.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 10-Oct-2014
 * \brief Implementation of the cylindrical constraint
 */ 

#include "constraint_cylinder.hpp"

/*! Force particle to be confined to the surface of a cylinder and
 *  its velocity to be tangent to it. We assume that during the integration 
 *  step particles do not move to much away from the surface and we project 
 *  them down to the surface of the cylinder
 *  \param p Particle which is to be projected onto the cylinder
 */
void ConstraintCylinder::enforce(Particle& p)
{
  bool apply = false;
  if (m_group == "all")
    apply = true;
  else
    apply = (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end());
  if (apply)
  {
    bool periodic = m_system->get_periodic();
    BoxPtr box = m_system->get_box();
    double x = p.x, y = p.y;
    double R = sqrt(x*x + y*y);
    double s = m_r/R;
    // Scale back to the surface
    p.x *= s; p.y *= s; 
    // Compute unit normal
    double Nx = p.x/m_r, Ny = p.y/m_r;
    // compute v.N
    double v_dot_N = p.vx*Nx + p.vy*Ny;
    // compute n.N
    double n_dot_N = p.nx*Nx + p.ny*Ny;
    // Project velocity onto tangent plane
    p.vx -= v_dot_N*Nx; p.vy -= v_dot_N*Ny;
    // Project director onto tangent plane
    p.nx -= n_dot_N*Nx; p.ny -= n_dot_N*Ny;
    // normalize director
    double inv_len = 1.0/sqrt(p.nx*p.nx + p.ny*p.ny + p.nz*p.nz);
    p.nx *= inv_len;  p.ny *= inv_len;  p.nz *= inv_len;
    if (periodic)
    {
      if (p.z > box->zhi) p.z -= box->Lz;
      else if (p.z < box->zlo) p.z += box->Lz;
    }
  }
}

/*! Rescale cylinder radius and make sure that all particles are still on it.
 *  Rescaling is done only at certain steps and only if rescale 
 *  factor is not equal to 1.
 */
bool ConstraintCylinder::rescale()
{
  if (m_rescale != 1.0)
  {
    int step = m_system->get_step();
    if (step % m_rescale_freq == 0 && step <= m_rescale_steps)
    {
      m_r *= m_scale;
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
