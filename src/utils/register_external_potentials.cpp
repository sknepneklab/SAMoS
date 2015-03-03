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
 * \file regster_external_potentials.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all external potentials with the class factory.
*/

#include "register.hpp"

void register_external_potentials(ExternalPotentialMap& external_potentials)
{
  // Register gravity to the external potential class factory
  external_potentials["gravity"] = boost::factory<ExternalGravityPotentialPtr>();
}