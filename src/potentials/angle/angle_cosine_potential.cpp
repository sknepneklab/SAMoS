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
 * \file angle_cosine_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com, modified by Amit Das, ncbs, India, dosamit@gmail.com
 * \date 05-Nov-2014
 * \brief Implementation of AngleCosinePotential class
 */ 

#include "angle_cosine_potential.hpp"

void AngleCosinePotential::compute()
{
  int Nangles = m_system->num_angles();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double k = m_k;
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < Nangles; i++)
  {
    Angle& a = m_system->get_angle(i);
    Particle& pi = m_system->get_particle(a.i);
    Particle& pj = m_system->get_particle(a.j);
    Particle& pk = m_system->get_particle(a.k);
    double dx1 = pi.x - pj.x, dy1 = pi.y - pj.y, dz1 = pi.z - pj.z;
    double dx2 = pk.x - pj.x, dy2 = pk.y - pj.y, dz2 = pk.z - pj.z;
    if (periodic)
    {
      if (dx1 > box->xhi) dx1 -= box->Lx;
      else if (dx1 < box->xlo) dx1 += box->Lx;
      if (dy1 > box->yhi) dy1 -= box->Ly;
      else if (dy1 < box->ylo) dy1 += box->Ly;
      if (dz1 > box->zhi) dz1 -= box->Lz;
      else if (dz1 < box->zlo) dz1 += box->Lz;
      if (dx2 > box->xhi) dx2 -= box->Lx;
      else if (dx2 < box->xlo) dx2 += box->Lx;
      if (dy2 > box->yhi) dy2 -= box->Ly;
      else if (dy2 < box->ylo) dy2 += box->Ly;
      if (dz2 > box->zhi) dz2 -= box->Lz;
      else if (dz2 < box->zlo) dz2 += box->Lz;
    }
    if (m_has_angle_params)
    {
      k = m_angle_params[a.type-1].k;
     }
      
    double r_sq_1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
    double r_sq_2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
    double r_1 = sqrt(r_sq_1);
    double r_2 = sqrt(r_sq_2);
    double c = dx1*dx2 + dy1*dy2 + dz1*dz2;
    c /= r_1*r_2;
    
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    
    m_potential_energy += k*(1.0+c);
    
    double aa = k;
    double a11 = aa*c / r_sq_1;
    double a12 = -aa / (r_1*r_2);
    double a22 = aa*c / r_sq_2;

    double fi_x = a11*dx1 + a12*dx2;
    double fi_y = a11*dy1 + a12*dy2;
    double fi_z = a11*dz1 + a12*dz2;
    
    double fk_x = a22*dx2 + a12*dx1;
    double fk_y = a22*dy2 + a12*dy1;
    double fk_z = a22*dz2 + a12*dz1;
    
    pi.fx += fi_x; pi.fy += fi_y; pi.fz += fi_z;
    
    pj.fx += -(fi_x+fk_x); pj.fy += -(fi_y+fk_y); pj.fz += -(fi_z+fk_z);
    
    pk.fx += fk_x; pk.fy += fk_y; pk.fz += fk_z;
    
  }
}
