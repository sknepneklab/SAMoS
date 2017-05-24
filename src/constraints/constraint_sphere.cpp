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
 * \file constraint_sphere.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 24-Oct-2013
 * \brief Implementation of the spherical constraint
 */ 

#include "constraint_sphere.hpp"

/*! Force particle to be confined to the surface of the sphere and
 *  its velocity to be tangent to it. We assume that during the integration 
 *  step particles do not move to much away from the surface and we project 
 *  them down using simple radial projection.
 *  \param p Particle which is to be projected onto the sphere
 */
void ConstraintSphere::enforce(Particle& p)
{
  bool apply = false;
  if (m_group == "all")
    apply = true;
  else
    apply = (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end());
  if (apply)
  {
    double x = p.x, y = p.y, z = p.z;
    double R = sqrt(x*x + y*y + z*z);
    double s = m_r/R;
    // Scale back to the surface
    p.x *= s; p.y *= s; p.z *= s;
    // Compute unit normal
    double Nx = p.x/m_r, Ny = p.y/m_r, Nz = p.z/m_r;
    // compute v.N
    double v_dot_N = p.vx*Nx + p.vy*Ny + p.vz*Nz;
    // compute n.N
    double n_dot_N = p.nx*Nx + p.ny*Ny + p.nz*Nz;
    // Project velocity onto tangent plane
    p.vx -= v_dot_N*Nx; p.vy -= v_dot_N*Ny; p.vz -= v_dot_N*Nz;
    // Project director onto tangent plane
    p.nx -= n_dot_N*Nx; p.ny -= n_dot_N*Ny; p.nz -= n_dot_N*Nz;
    // normalize director
    double inv_len = 1.0/sqrt(p.nx*p.nx + p.ny*p.ny + p.nz*p.nz);
    p.nx *= inv_len;  p.ny *= inv_len;  p.nz *= inv_len;
    // Set particle normal
    p.Nx = Nx; p.Ny = Ny; p.Nz = Nz;
    // Project all forces onto tangent plane
    if (m_system->record_force_type())
    {
      map<string,ForceType>& force_type = p.get_force_type();
      for (map<string,ForceType>::iterator it = force_type.begin(); it != force_type.end(); it++)
      {
        double fx = (*it).second.fx, fy = (*it).second.fy, fz = (*it).second.fz; 
        double f_dot_N = fx*Nx + fy*Ny + fz*Nz;
        (*it).second.fx -= f_dot_N*Nx; 
        (*it).second.fy -= f_dot_N*Ny;
        (*it).second.fz -= f_dot_N*Nz;
      }
    }
  }
}

/*! Rescale sphere size and make sure that all particles are still on it.
 *  Rescaling is done only at certain steps and only if rescale 
 *  factor is not equal to 1.
 */
bool ConstraintSphere::rescale()
{
  if (m_rescale != 1.0)
  {
    int step = m_system->get_step();
    if ((step % m_rescale_freq == 0) && (step < m_rescale_steps))
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

