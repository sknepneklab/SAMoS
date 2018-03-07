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
 * \file external_harmonic_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-Aug-2015
 * \brief Implementation of the external harmonic potential
 */ 

#include "external_harmonic_potential.hpp"

/*! Apply external harmonic potential to all particles */
void ExternalHarmonicPotential::compute()
{
  int N = m_system->size();
  double k = m_k;
  double z0 = m_z0;
  
  m_potential_energy = 0.0;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    if (m_has_params)
    {
      k = m_type_params[p.get_type()-1]["k"];
      z0 = m_type_params[p.get_type()-1]["z0"];
    }
    double dz = p.z - z0;
    m_potential_energy += 0.5*k*dz*dz;
    p.fz -= k*dz;
  }
}
