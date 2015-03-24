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
 * \file regster_pair_potentials.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all pair potentials with the class factory.
*/

#include "register.hpp"

void register_pair_potentials(PairPotentialMap& pair_potentials)
{
  // Register Lennard-Jones pair potential with the pair potentials class factory
  pair_potentials["lj"] = boost::factory<PairLJPotentialPtr>();
  // Register Coulomb pair potential with the pair potentials class factory
  pair_potentials["coulomb"] = boost::factory<PairCoulombPotentialPtr>();
  // Register soft pair potential with the pair potentials class factory
  pair_potentials["soft"] = boost::factory<PairSoftPotentialPtr>();
  // Register Gaussian pair potential with the pair potentials class factory
  pair_potentials["gaussian"] = boost::factory<PairGaussianPotentialPtr>();  
  // Register Morse pair potential with the pair potentials class factory
  pair_potentials["morse"] = boost::factory<PairMorsePotentialPtr>();  
  // Register active pair potential with the pair potentials class factory
  pair_potentials["active"] = boost::factory<PairActivePotentialPtr>(); 
  // Register soft rod pair potential with the pair potentials class factory
  pair_potentials["rod"] = boost::factory<PairRodPotentialPtr>(); 
}
