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
 * \file constraint_slab.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Sep-2015
 * \brief Implementation of the slab constraint
 */ 

#include "constraint_slab.hpp"

/*! Force all particles to be confined move between planes at \f$ z = z_{lo} \f$ and 
 *  \f$ z = z_{hi} \f$.
 */
void ConstraintSlab::enforce(Particle& p)
{
  bool apply = false;
  if (m_group == "all")
    apply = true;
  else
    apply = (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end());
  if (apply)
  {
    if (p.z < m_z_lo)
    {
      p.z = m_z_lo;
      p.vz = -p.vz;
      p.nz = -p.nz;
    }
    if (p.z > m_z_hi)
    {
      p.z = m_z_hi;
      p.vz = -p.vz;
      p.nz = -p.nz;
    }
  }
  if (m_system->get_periodic())
    m_system->enforce_periodic(p);
}



