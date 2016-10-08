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
 * \file regster.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Auxiliary functions for registering different simulation 
 *        quantities with the corresponding class factory
*/

#ifndef __REGISTER_HPP__
#define __REGISTER_HPP__

#include "factory_types.hpp"
#include "samos.hpp"

// Register constraints 
void register_constraints(ConstraintMap&);
//Register pair potentials
void register_pair_potentials(PairPotentialMap&);
//Register external potentials
void register_external_potentials(ExternalPotentialMap&);
//Register integrators
void register_integrators(IntegratorMap&);
//Register pair aligners
void register_pair_aligners(PairAlignerMap&);
//Register external aligners
void register_external_aligners(ExternalAlignerMap&);
//Register populations
void register_populations(PopulationMap&);
//Register bond potentials
void register_bond_potentials(BondPotentialMap&);
//Register angle potentials
void register_angle_potentials(AnglePotentialMap&);
//Register values
void register_values(ValueMap&);

#endif
