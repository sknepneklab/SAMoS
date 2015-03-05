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
 * \file regster_values.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all values with the class factory.
*/

#include "register.hpp"

void register_values(ValueMap& values)
{
  // Register constant value control the class factory
  values["constant"] = boost::factory<ValueConstantPtr>();
  // Register linear value control the class factory
  values["linear"] = boost::factory<ValueLinearPtr>();
}