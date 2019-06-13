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

#include "register.hpp"

void register_bond_potentials(BondPotentialMap& bond_potentials)
{
  // Register harmonic bond potential with the class factory
  bond_potentials["harmonic"] = factory<BondHarmonicPotentialPtr>();
  // Register active bond force with the class factory
  bond_potentials["active"] = factory<BondActiveForcePtr>();
  // Register FENE bond force with the class factory
  bond_potentials["fene"] = factory<BondHarmonicPotentialPtr>();
}
