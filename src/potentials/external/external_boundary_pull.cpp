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
 * \file external_boundary_pull.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Sep-2016
 * \brief Implementation of the self propulsion for active particles.
 */ 

#include "external_boundary_pull.hpp"

/*! Apply pulling force onto each boundary particle. */
void ExternalBoundaryPull::compute()
{
  int N = m_system->size();
  double alpha = m_alpha;

  // First we compute centre of mass to be albe to determine outward direction
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
      double force_sign = 1.0;
      double X = pi.x - xcm,  Y = pi.y - ycm,  Z = pi.z - zcm;
      double len_R = sqrt(X*X + Y*Y + Z*Z);
      X /= len_R;  Y /= len_R;  Y /= len_R;
      Particle& pj = m_system->get_particle(pi.boundary_neigh[0]);
      Particle& pk = m_system->get_particle(pi.boundary_neigh[1]);
      double xji = pj.x - pi.x, yji = pj.y - pi.y, zji = pj.z - pi.z;
      double len_ji = sqrt(xji*xji + yji*yji + zji*zji);
      xji /= len_ji;  yji /= len_ji;  zji /= len_ji; 
      double xki = pk.x - pi.x, yki = pk.y - pi.y, zki = pk.z - pi.z;
      double len_ki = sqrt(xki*xki + yki*yki + zki*zki);
      xki /= len_ki;  yki /= len_ki;  zki /= len_ki;
      double x = -(xji+xki), y = -(yji+yki), z = -(zji+zki);
      double len_r = sqrt(x*x + y*y + z*z);
      x /= len_r;  y /= len_r;  z /= len_r;
      // compute dot product with radius the vector connecting center of mass and pi
      if ((x*X + y*Y + z*Z) < 0.0)
        force_sign = -1.0;
      double factor = force_sign*alpha;
      pi.fx += factor*x;
      pi.fy += factor*y;
      pi.fz += factor*z;
    }
  }
}
