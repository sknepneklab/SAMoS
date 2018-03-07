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
 * \file external_ajnematic_aligner.cpp
 * \author Silke Henkes, silkehenkes@gmail.com
 * \date 16-Oct-2014
 * \brief Declaration of ExternalAJNematicAlign class
 */ 


#include "external_ajnematic_aligner.hpp"

void ExternalAJNematicAlign::compute()
{
  int N = m_system->size();
  double tau = m_tau;
  
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    // The active jamming torque for nematic systems is modelled as:
    // 1/tau sin(2(\theta_v-theta_i))
    // or in vectorial: -2/tau (n_i.v_i)(n_i \times v_i)
    // note: factor 2 comes from the expansion sin(2x) = 2sin(x)cos(x)
    double ni_dot_vi = pi.nx*pi.vx + pi.ny*pi.vy + pi.nz*pi.vz;
    double tau_x = pi.ny*pi.vz - pi.nz*pi.vy;
    double tau_y = pi.nz*pi.vx - pi.nx*pi.vz;
    double tau_z = pi.nx*pi.vy - pi.ny*pi.vx;
    
     // Do normalisation here if useful
    if (m_normalise)
    {
      double vnorm = sqrt(pi.vx*pi.vx + pi.vy*pi.vy + pi.vz*pi.vz);
      if (vnorm > 0.0)
      {
        tau_x /= vnorm;
	      tau_y /= vnorm;
	      tau_z /= vnorm;
      }
    }
    
    if (m_has_params)
      tau = m_type_params[pi.get_type()-1].tau;
    
    pi.tau_x += ni_dot_vi*tau_x/tau;
    pi.tau_y += ni_dot_vi*tau_y/tau;
    pi.tau_z += ni_dot_vi*tau_z/tau;
  }
}
