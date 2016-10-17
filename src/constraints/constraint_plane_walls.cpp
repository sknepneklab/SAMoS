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
 * \file constraint_plane_walls.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 13-Aug-2014
 * \brief Implementation of the planar walls constraint
 */ 

#include "constraint_plane_walls.hpp"

/*! Force all particles to be confined to the surface of the xy plane 
 *  by setting z = 0, vz = 0, fz = 0.
 *  At the walls, particle is stopped at the wall by setting \f$ v_x = 0 \f$
 *  and x-coordinate is set back to l or -l.
 *  \param p particle to project onto the sphere
 */
void ConstraintPlaneWalls::enforce(Particle& p)
{
  bool apply = false;
  if (m_group == "all")
    apply = true;
  else
    apply = (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end());
  if (apply)
  {
    bool periodic = m_system->get_periodic();
    double ylo = m_system->get_box()->ylo, yhi = m_system->get_box()->yhi;
    double ly = yhi - ylo;
    p.z = 0.0;
    p.vz = 0.0;
    p.fz = 0.0;
    // Check periodic boundary conditions 
    if (periodic)
    {
      if (p.y <= ylo) p.y += ly;
      else if (p.y >= yhi) p.y -= ly;
    }
    if (p.x >= m_l)
    {
      p.x = m_l;
      p.vx = -p.vx;
    }
    else if (p.x <= -m_l)
    {
      p.x = -m_l;
      p.vx = -p.vx;
    }
  }
}

/*! Rotate director of a particle around the normal vector (z axis)
 *  \note This function assumes that the particle has already been
 *  projected onto the plane and that its director is also in plane
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintPlaneWalls::rotate_director(Particle& p, double phi)
{
  bool apply = false;
  if (m_group == "all")
    apply = true;
  else
    apply = (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end());
  if (apply)
  {
    // Sine and cosine of the rotation angle
    double c = cos(phi), s = sin(phi);
    // Rotation matrix around z axis
    double Rxx = c, Rxy = -s;
    double Ryx = s, Ryy = c;
    // Apply rotation matrix
    double nx = Rxx*p.nx + Rxy*p.ny;
    double ny = Ryx*p.nx + Ryy*p.ny;
    double len = sqrt(nx*nx + ny*ny);
    // Update particle director
    p.nx = nx/len;
    p.ny = ny/len;
  }
}

/*! Rotate velocity of a particle around the normal vector (z axis)
 *  \note This function assumes that the particle has already been
 *  projected onto the plane and that its director is also in plane
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintPlaneWalls::rotate_velocity(Particle& p, double phi)
{
  bool apply = false;
  if (m_group == "all")
    apply = true;
  else
    apply = (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end());
  if (apply)
  {
    // Sine and cosine of the rotation angle
    double c = cos(phi), s = sin(phi);
    // Rotation matrix around z axis
    double Rxx = c, Rxy = -s;
    double Ryx = s, Ryy = c;
    // Apply rotation matrix
    double vx = Rxx*p.vx + Rxy*p.vy;
    double vy = Ryx*p.vx + Ryy*p.vy;
    // Update particle director
    p.vx = vx;
    p.vy = vy;
  }
}

/*! Project particle torque onto the normal vector (z axis). The assumption here is that 
 *  the particle's director and velocity are all already in the xy plane
 *  and that it is constrained to the plane.
 *  \param p Particle whose torque to project
*/ 
double ConstraintPlaneWalls::project_torque(Particle& p)
{
  bool apply = false;
  if (m_group == "all")
    apply = true;
  else
    apply = (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end());
  
  if (apply)
    return p.tau_z;  
  else
    return 0.0;
}
