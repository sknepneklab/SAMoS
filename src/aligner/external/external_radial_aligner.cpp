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
 * \file external_radial_aligner.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Feb-2018
 * \brief Declaration of ExternalRadialAlign class
 */ 


#include "external_radial_aligner.hpp"

void ExternalRadialAlign::compute()
{
  int N = m_system->size();
  double J = m_J;
  double pos_x = m_pos_x, pos_y = m_pos_y;
  double x_c = 0.0, y_c = 0.0;
  if (m_use_cm)
  {
    int cnt = 0;
    for (int i = 0; i < N; i++)
    {
      Particle& pi = m_system->get_particle(i);
      if (pi.in_tissue)
      {
        x_c += pi.x;
        y_c += pi.y;
        cnt++;
      }
    }
    x_c /= cnt;  y_c /= cnt;
  } 

  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (pi.in_tissue)
    {
      if (m_has_params)
      {
        J = m_type_params[pi.get_type()-1].J;
        pos_x = m_type_params[pi.get_type()-1].pos_x;
        pos_y = m_type_params[pi.get_type()-1].pos_y;
      }
      double X = pos_x + x_c, Y = pos_y + y_c;
      double x = pi.x - X, y = pi.y - Y;
      if ((x != 0.0) || (y != 0.0))
      {
        double r = sqrt(x*x + y*y);
        double fact = -J/r;
        double tau_z = x*pi.ny - y*pi.nx;
        pi.tau_z += fact*tau_z;
      }
    }
  }
}
