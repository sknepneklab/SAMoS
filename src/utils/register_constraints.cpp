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
 * \file regster_constraints.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all constraints with the class factory.
*/

#include "register.hpp"

void register_constraints(ConstraintMap& constraints)
{
  // Register spherical constraint with the constraint class factory
  constraints["sphere"] = boost::factory<ConstraintSpherePtr>();
  // Register xy plane constraint with the constraint class factory
  constraints["plane"] = boost::factory<ConstraintPlanePtr>();
  // Register plane walls constraint with the constraint class factory
  constraints["walls"] = boost::factory<ConstraintPlaneWallsPtr>();
  // Register cylindrical constraint with the constraint class factory
  constraints["cylinder"] = boost::factory<ConstraintCylinderPtr>();
  // Register peanut constraint with the constraint class factory
  constraints["peanut"] = boost::factory<ConstraintPeanutPtr>();
  // Register torus constraint with the constraint class factory
  constraints["torus"] = boost::factory<ConstraintTorusPtr>();
  // Register ellipsoid constraint with the constraint class factory
  constraints["ellipsoid"] = boost::factory<ConstraintEllipsoidPtr>();
  // Register gyroid constraint with the constraint class factory
  constraints["gyroid"] = boost::factory<ConstraintGyroidPtr>();
  // Register actomyo constraint with the constraint class factory
  constraints["actomyo"] = boost::factory<ConstraintActomyoPtr>();
  // Register hourglass constraint with the constraint class factory
  constraints["hourglass"] = boost::factory<ConstraintHourglassPtr>();
  // Register Gaussian bump constraint with the constraint class factory
  constraints["gaussian_bump"] = boost::factory<ConstraintGaussianBumpPtr>();
  // Register dummy constraint with the constraint class factory
  constraints["none"] = boost::factory<ConstraintNonePtr>();
}
