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
  constraints["sphere"] = factory<ConstraintSpherePtr>();
  // Register xy plane constraint with the constraint class factory
  constraints["plane"] = factory<ConstraintPlanePtr>();
  // Register plane walls constraint with the constraint class factory
  constraints["walls"] = factory<ConstraintPlaneWallsPtr>();
  // Register cylindrical constraint with the constraint class factory
  constraints["cylinder"] = factory<ConstraintCylinderPtr>();
  // Register peanut constraint with the constraint class factory
  constraints["peanut"] = factory<ConstraintPeanutPtr>();
  // Register torus constraint with the constraint class factory
  constraints["torus"] = factory<ConstraintTorusPtr>();
  // Register ellipsoid constraint with the constraint class factory
  constraints["ellipsoid"] = factory<ConstraintEllipsoidPtr>();
  // Register gyroid constraint with the constraint class factory
  constraints["gyroid"] = factory<ConstraintGyroidPtr>();
  // Register actomyo constraint with the constraint class factory
  constraints["actomyo"] = factory<ConstraintActomyoPtr>();
  // Register hourglass constraint with the constraint class factory
  constraints["hourglass"] = factory<ConstraintHourglassPtr>();
  // Register Gaussian bump constraint with the constraint class factory
  constraints["gaussian_bump"] = factory<ConstraintGaussianBumpPtr>();
  // Register dummy constraint with the constraint class factory
  constraints["none"] = factory<ConstraintNonePtr>();
  // Register constraint on the surface with tetraherdal symmetry with the constraint class factory
  constraints["tetrahedron"] = factory<ConstraintTetrahedronPtr>();
  // Register constraint to move in a slab between to planes parallel to xy plane with the constraint class factory
  constraints["slab"] = factory<ConstraintSlabPtr>();
}
