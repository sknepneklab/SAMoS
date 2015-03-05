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
 * \file regster_integrators.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all integrators with the class factory.
*/

#include "register.hpp"

void register_integrators(IntegratorMap& integrators)
{
  // Register Brownian dynamics integrator with the integrators class factory
  integrators["brownian"] = boost::factory<IntegratorBrownianPtr>();
  // Register Vicsek dynamics integrator with the integrators class factory
  integrators["vicsek"] = boost::factory<IntegratorVicsekPtr>();
  // Register NVE integrator with the integrators class factory
  integrators["nve"] = boost::factory<IntegratorNVEPtr>();
  // Register nematic integrator with the integrators class factory
  integrators["nematic"] = boost::factory<IntegratorNematicPtr>();
  // Register actomyo integrator with the integrators class factory
  integrators["actomyo"] = boost::factory<IntegratorActomyoPtr>();
}