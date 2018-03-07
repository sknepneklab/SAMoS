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


#include "external_kenotaxis_aligner.hpp"

void ExternalKenotaxisAlign::compute()
{
  int N = m_system->size();
  double J = m_J;
  Mesh& mesh = m_system->get_mesh();


  // First we compute centre of mass to be able to determine outward direction
  double xcm = 0.0, ycm = 0.0, zcm = 0.0;
  int cnt = 0;
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (pi.in_tissue)
    {
      xcm += pi.x;
      ycm += pi.y;
      zcm += pi.z;
      cnt++;
    }
  }
  xcm /= cnt;  ycm /= cnt;  zcm /= cnt;
  
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (pi.boundary)
    {
      // Compute the direction perpendicular to the boundary.
      // We cheat by using the centre of mass position to deterimine the outside and inside of the tissue.

      // vector (X, Y, Z) points from centre of mass to particle i
      double force_sign = 1.0;
      double X = pi.x - xcm,  Y = pi.y - ycm,  Z = pi.z - zcm;
      double len_R = sqrt(X*X + Y*Y + Z*Z);
      X /= len_R;  Y /= len_R;  Z /= len_R;
      // get relative vectors of neighbouring boundary particles
      Particle& pj = m_system->get_particle(pi.boundary_neigh[0]);
      Particle& pk = m_system->get_particle(pi.boundary_neigh[1]);
      double xji = pj.x - pi.x, yji = pj.y - pi.y, zji = pj.z - pi.z;
      double len_ji = sqrt(xji*xji + yji*yji + zji*zji);
      xji /= len_ji;  yji /= len_ji;  zji /= len_ji; 
      double xki = pk.x - pi.x, yki = pk.y - pi.y, zki = pk.z - pi.z;
      double len_ki = sqrt(xki*xki + yki*yki + zki*zki);
      xki /= len_ki;  yki /= len_ki;  zki /= len_ki;
      // compute vector perpendicular to r_kj
      double x = -(yji-yki), y = (xji-xki), z = 0.0; // give up with three dimensions, z=0
      double len_r = sqrt(x*x + y*y + z*z);
      x /= len_r;  y /= len_r;  z /= len_r;
      // compute dot product with vector connecting center of mass and (x,y,z)
      if ((x*X + y*Y + z*Z) < 0.0)
        force_sign = -1.0;
      double J_factor = force_sign*J;
      double tau_z = pi.nx*y - pi.ny*x;
      //if (m_has_params)
        //J = m_type_params[pi.get_type()-1].J
        //
      
      pi.tau_x += 0.0;
      pi.tau_y += 0.0;
      pi.tau_z += J_factor*tau_z;
    }
  }
}
