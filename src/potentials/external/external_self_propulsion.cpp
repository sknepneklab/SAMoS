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
 * \file external_self_propulsion.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Dec-2015
 * \brief Implementation of the self propulsion for active particles.
 */ 

#include "external_self_propulsion.hpp"

/*! Apply self propulsion to all particles */
void ExternalSelfPropulsion::compute()
{
  int N = m_system->size();
  double alpha = m_alpha;
  
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    if (!(m_exclude_boundary && p.boundary))
    {
      if (m_has_params)
        alpha = m_type_params[p.get_type()-1]["alpha"];
      p.fx += alpha*p.nx;
      p.fy += alpha*p.ny;
      p.fz += alpha*p.nz;
    }
  }
}
