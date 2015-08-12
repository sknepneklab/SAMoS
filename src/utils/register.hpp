/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

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
