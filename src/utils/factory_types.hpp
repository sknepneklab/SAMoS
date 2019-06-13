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
 * \file factory_types.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Define class factory data types
*/

#ifndef __FACTORY_TYPES_HPP__
#define __FACTORY_TYPES_HPP__

#include <map>
#include <functional>
#include <string>


#include "samos.hpp"
#include "class_factory.hpp"


typedef std::function< ConstraintPtr(SystemPtr, MessengerPtr, pairs_type&) > constraint_factory;                                                                        //!< Type that handles all constraints 
typedef std::function< PairPotentialPtr(SystemPtr, MessengerPtr, NeighbourListPtr, ValuePtr, pairs_type&) > pair_potential_factory;                                     //!< Type that handles all pair potentials
typedef std::function< ExternalPotentialPtr(SystemPtr, MessengerPtr, pairs_type&) > external_potential_factory;                                                         //!< Type that handles all external potentials
typedef std::function< IntegratorPtr(SystemPtr, MessengerPtr, PotentialPtr, AlignerPtr, NeighbourListPtr, ConstrainerPtr, ValuePtr, pairs_type&) > integrator_factory;  //!< Type that handles all integrators
typedef std::function< PairAlignPtr(SystemPtr, MessengerPtr, NeighbourListPtr, pairs_type&) > pair_aligner_factory;                                                     //!< Type that handles all pairwise alignment
typedef std::function< ExternalAlignPtr(SystemPtr, MessengerPtr, pairs_type&) > external_aligner_factory;                                                               //!< Type that handles all external alignment
typedef std::function< PopulationPtr(SystemPtr, MessengerPtr, pairs_type&) > population_factory;                                                                        //!< Type that handles all populations
typedef std::function< BondPotentialPtr(SystemPtr, MessengerPtr, pairs_type&) > bond_potential_factory;                                                                 //!< Type that handles all bond potentials
typedef std::function< AnglePotentialPtr(SystemPtr, MessengerPtr, pairs_type&) > angle_potential_factory;                                                               //!< Type that handles all angle potentials
typedef std::function< ValuePtr(MessengerPtr, pairs_type&) > value_factory;                                                                                             //!< Type that handles all value control objects


typedef std::map<std::string, constraint_factory> ConstraintMap;
typedef std::map<std::string, pair_potential_factory> PairPotentialMap;
typedef std::map<std::string, external_potential_factory> ExternalPotentialMap;
typedef std::map<std::string, integrator_factory> IntegratorMap;
typedef std::map<std::string, pair_aligner_factory> PairAlignerMap;
typedef std::map<std::string, external_aligner_factory> ExternalAlignerMap;
typedef std::map<std::string, population_factory> PopulationMap;
typedef std::map<std::string, bond_potential_factory> BondPotentialMap;
typedef std::map<std::string, angle_potential_factory> AnglePotentialMap;
typedef std::map<std::string, value_factory> ValueMap;


#endif
