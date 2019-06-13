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
 * \file regster_external_potentials.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all external potentials with the class factory.
*/

#include "register.hpp"

void register_external_potentials(ExternalPotentialMap& external_potentials)
{
  // Register gravity to the external potential class factory
  external_potentials["gravity"] = factory<ExternalGravityPotentialPtr>();
  // Register harmonic potential to the external potential class factory
  external_potentials["harmonic"] = factory<ExternalHarmonicPotentialPtr>();
  // Register self propulsion to the external potential class factory
  external_potentials["self_propulsion"] = factory<ExternalSelfPropulsionPtr>();
  // Register boundary pull to the external potential class factory
  external_potentials["boundary_pull"] = factory<ExternalBoundaryPullPtr>();
}
