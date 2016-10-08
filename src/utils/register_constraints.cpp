/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

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
  // Register constraint on the surface with tetraherdal symmetry with the constraint class factory
  constraints["tetrahedron"] = boost::factory<ConstraintTetrahedronPtr>();
  // Register constraint to move in a slab between to planes parallel to xy plane with the constraint class factory
  constraints["slab"] = boost::factory<ConstraintSlabPtr>();
}
