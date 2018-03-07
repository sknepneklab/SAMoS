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
 * \file external_gravity_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 23-Oct-2013
 * \brief Implementation of the external gravitational potential
 */

#include "external_gravity_potential.hpp"

/*! Apply external gravitational potential to all particles */
void ExternalGravityPotential::compute()
{
  int N = m_system->size();
  double g = m_g;

  m_potential_energy = 0.0;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    if (m_has_params)
      g = m_type_params[p.get_type()-1]["g"];
    m_potential_energy += g*p.z;
    p.fz -= g;
  }
}
