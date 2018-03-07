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
 * \file external_tangent_aligner.cpp
 * \author Silke Henkes, silkehenkes@gmail.com
 * \date 16-Oct-2014
 * \brief Declaration of ExternaltangentAlign class
 */ 


#include "external_tangent_aligner.hpp"

void ExternalTangentAlign::compute()
{
  int N = m_system->size();
  double tau = m_tau;
  int Nbonds = m_system->num_bonds();
    
  for  (int i = 0; i < Nbonds; i++)
  {
    Bond& b = m_system->get_bond(i);
    Particle& pi = m_system->get_particle(b.i);
    Particle& pj = m_system->get_particle(b.j);
    double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;  // tangent along bond between i and j
    m_system->apply_periodic(dx,dy,dz);
    double l = sqrt(dx*dx + dy*dy + dz*dz);     
    dx /= l; dy /= l; dz /= l;              // normalise tangent vector
    // Compute torque between director of particle i and tangent
    double tau_x = pi.ny*dz - pi.nz*dy;
    double tau_y = pi.nz*dx - pi.nx*dz;
    double tau_z = pi.nx*dy - pi.ny*dx;
    if (m_has_params)
      tau = m_type_params[pi.get_type()-1].tau;
    pi.tau_x += tau_x/tau;
    pi.tau_y += tau_y/tau;
    pi.tau_z += tau_z/tau;
    // Compute torque between director of particle j and tangent
    tau_x = pj.ny*dz - pj.nz*dy;
    tau_y = pj.nz*dx - pj.nx*dz;
    tau_z = pj.nx*dy - pj.ny*dx;
    if (m_has_params)
      tau = m_type_params[pj.get_type()-1].tau;
    pj.tau_x += tau_x/tau;
    pj.tau_y += tau_y/tau;
    pj.tau_z += tau_z/tau;
  }
}
