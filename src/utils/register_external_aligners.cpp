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
 * \file regster_external_aligners.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all external aligners with the class factory.
*/

#include "register.hpp"

void register_external_aligners(ExternalAlignerMap& external_aligners)
{
  // Register active jamming polar external alignment to the external align class factory
  external_aligners["aj"] = boost::factory<ExternalAJPolarAlignPtr>();
  // Register external alignment a vector filed to the external align class factory
  external_aligners["field"] = boost::factory<ExternalFieldAlignPtr>();
}
