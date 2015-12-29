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
  // Register Brownian rod dynamics integrator with the integrators class factory
  integrators["brownian_rod"] = boost::factory<IntegratorBrownianRodPtr>();
  // Register Brownian dynamics integrator for particle position with the integrators class factory
  integrators["brownian_pos"] = boost::factory<IntegratorBrownianPosPtr>();
  // Register Brownian dynamics integrator for alignment with the integrators class factory
  integrators["brownian_align"] = boost::factory<IntegratorBrownianAlignPtr>();
  // Register Langevin dynamics integrator for particle positions with the integrators class factory
  integrators["langevin"] = boost::factory<IntegratorLangevinPtr>();
}