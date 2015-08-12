/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014, 2015
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file regster_bond_potentials.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all bond potentials with the class factory.
*/

#include "register.hpp"

void register_bond_potentials(BondPotentialMap& bond_potentials)
{
  // Register harmonic bond potential with the class factory
  bond_potentials["harmonic"] = boost::factory<BondHarmonicPotentialPtr>();
  // Register active bond force with the class factory
  bond_potentials["active"] = boost::factory<BondActiveForcePtr>();
}