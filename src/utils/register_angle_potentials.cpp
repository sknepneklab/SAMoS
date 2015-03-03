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
 * \file regster_angle_potentials.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all angle potentials with the class factory.
*/

#include "register.hpp"

void register_angle_potentials(AnglePotentialMap& angle_potentials)
{
  // Register harmonic angle potential with the class factory
  angle_potentials["harmonic"] = boost::factory<AngleHarmonicPotentialPtr>();
    // Register cosine angle potential with the class factory
  angle_potentials["cosine"] = boost::factory<AngleCosinePotentialPtr>();
}