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
 * \file regster_populations.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all populations with the class factory.
*/

#include "register.hpp"

void register_populations(PopulationMap& populations)
{
  // Register random population control with the class factory
  populations["random"] = boost::factory<PopulationRandomPtr>();
  // Register density population control with the class factory
  populations["density"] = boost::factory<PopulationDensityPtr>();
}