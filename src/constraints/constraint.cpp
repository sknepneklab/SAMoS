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
 * \file constraint.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 09-Feb-2015
 * \brief Implementation of the generic constraint
 */ 

#include "constraint.hpp"


/*! Force particle to be confined to the constraining surface and
 *  its velocity to be tangent to it. We assume that during the integration 
 *  step particles do not move to much away from the surface and we project 
 *  them down to the surface. This is a generic iterative method. Some constraints, 
 *  like sphere override this method to make it faster.
 *  \param p Particle which is to be projected onto the constraint
 */
void Constraint::enforce(Particle& p)
{
  bool apply = false;
  if (m_group == "all")
    apply = true;
  else
    apply = (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end());
  if (apply)
  {
    // Compute reference gradient (SHAKE method)
    double ref_gx, ref_gy, ref_gz;
    this->compute_gradient(p,ref_gx, ref_gy, ref_gz);
    
    bool satisfied;
    int iter = 0;
    while (iter++ < m_max_iter)
    {
      satisfied = true;
      // Compute unit gradient of the implicit function
      double gx, gy, gz;
      this->compute_gradient(p,gx,gy,gz);
      double s = gx*ref_gx + gy*ref_gy + gz*ref_gz;
      double g = this->constraint_value(p);
      double lambda = g/s;
      if (fabs(g) > m_tol)
        satisfied = false;
      if (satisfied) break;
      p.x -= lambda*ref_gx;
      p.y -= lambda*ref_gy;
      p.z -= lambda*ref_gz;
    }
      
    double Nx, Ny, Nz;
    this->compute_normal(p,Nx,Ny,Nz);
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
    m_system->enforce_periodic(p);
  }
}

/*! Rotate director of a particle around the normal vector
 *  \note This function assumes that the particle has already been
 *  projected onto the surface and that its director is laying in 
 *  the tangent plane.
 *  \param p Particle whose director to rotate
 *  \param phi angle by which to rotate it
*/
void Constraint::rotate_director(Particle& p, double phi)
{
  bool apply = false;
  if (m_group == "all")
    apply = true;
  else
    apply = (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end());
  if (apply)
  {
    double U, V, W;
    this->compute_normal(p,U,V,W);
    // Compute angle sins and cosines
    double c = cos(phi), s = sin(phi);
    // Compute new velocity coordinates
    double nx = U*(U*p.nx+V*p.ny+W*p.nz)*(1.0-c) + p.nx*c + (-W*p.ny+V*p.nz)*s;
    double ny = V*(U*p.nx+V*p.ny+W*p.nz)*(1.0-c) + p.ny*c + ( W*p.nx-U*p.nz)*s;
    double nz = W*(U*p.nx+V*p.ny+W*p.nz)*(1.0-c) + p.nz*c + (-V*p.nx+U*p.ny)*s;
    double len = sqrt(nx*nx + ny*ny + nz*nz);
    // Update particle director (normalize it along the way to collect for any numerical drift that may have occurred)
    p.nx = nx/len;
    p.ny = ny/len;
    p.nz = nz/len;  
  }
}

/*! Rotate velocity of a particle around the normal vector
 *  \note This function assumes that the particle has already been
 *  projected onto the surface and that its velocity is laying in 
 *  the tangent plane.
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void Constraint::rotate_velocity(Particle& p, double phi)
{
  bool apply = false;
  if (m_group == "all")
    apply = true;
  else
    apply = (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end());
  if (apply)
  {
    double U, V, W;
    this->compute_normal(p,U,V,W);
    // Compute angle sins and cosines
    double c = cos(phi), s = sin(phi);
    // Compute new velocity coordinates
    double vx = U*(U*p.vx+V*p.vy+W*p.vz)*(1.0-c) + p.vx*c + (-W*p.vy+V*p.vz)*s;
    double vy = V*(U*p.vx+V*p.vy+W*p.vz)*(1.0-c) + p.vy*c + ( W*p.vx-U*p.vz)*s;
    double vz = W*(U*p.vx+V*p.vy+W*p.vz)*(1.0-c) + p.vz*c + (-V*p.vx+U*p.vy)*s;
    // Update particle velocity
    p.vx = vx;
    p.vy = vy;
    p.vz = vz;  
  }
}


/*! Project particle torque onto the normal vector. The assumption here is that 
 *  the particle's director and velocity are all already in the tangent plane 
 *  and that it is constrained to the surface.
 *  \param p Particle whose torque to project
*/ 
double Constraint::project_torque(Particle& p)
{
  bool apply = false;
  if (m_group == "all")
    apply = true;
  else
    apply = (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end());
  if (apply)
  {
    double Nx, Ny, Nz;
    this->compute_normal(p,Nx,Ny,Nz);
    return (p.tau_x*Nx + p.tau_y*Ny + p.tau_z*Nz);  
  }
  else 
    return 0.0;
}
