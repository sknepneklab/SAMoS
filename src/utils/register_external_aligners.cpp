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
 * \file regster_external_aligners.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all external aligners with the class factory.
*/

#include "register.hpp"

void register_external_aligners(ExternalAlignerMap& external_aligners)
{
  // Register active jamming polar external alignment to the external align class factory
  external_aligners["ajpolar"] = boost::factory<ExternalAJPolarAlignPtr>();
  // Register external alignment a vector filed to the external align class factory
  external_aligners["field"] = boost::factory<ExternalFieldAlignPtr>();
  // Register active jamming nematic external alignment to the external align class factory
  external_aligners["ajnematic"] = boost::factory<ExternalAJNematicAlignPtr>();
  // Register cell shape external alignment to the external align class factory
  external_aligners["cell_shape"] = boost::factory<ExternalShapeAlignPtr>();
}
