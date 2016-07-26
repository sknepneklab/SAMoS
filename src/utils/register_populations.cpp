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
  // Register growth population control with the class factory
  populations["grow"] = boost::factory<PopulationGrowthPtr>();
  // Register elongation population control with the class factory
  populations["elongate"] = boost::factory<PopulationElongationPtr>();
  // Register cell population control with the class factory
  populations["cell"] = boost::factory<PopulationCellPtr>();
  // Register actomyosin population control with the class factory
  populations["actomyosin"] = boost::factory<PopulationActomyosinPtr>();
  // Register actomyosin poisson population control with the class factory
  populations["actomyosin_poisson"] = boost::factory<PopulationActomyosinPoissonPtr>();
}
