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
  // Register Lennard-Jones rod pair potential with the pair potentials class factory
  pair_potentials["ljrod"] = boost::factory<PairLJRodPotentialPtr>(); 
  // Register soft attractive pair potential with the pair potentials class factory
  pair_potentials["soft_attractive"] = boost::factory<PairSoftAttractivePotentialPtr>();
  // Register vertex particle pair potential for tissues with the pair potentials class factory
  pair_potentials["vp"] = boost::factory<PairVertexParticlePotentialPtr>();
  // Register line tension pair potential for tissues with the pair potentials class factory
  pair_potentials["line_tension"] = boost::factory<PairLineTensionPotentialPtr>();
}
