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
 * \file regster_external_aligners.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all external aligners with the class factory.
*/

#include "register.hpp"

void register_external_aligners(ExternalAlignerMap& external_aligners)
{
  // Register active jamming polar external alignment to the external align class factory
  external_aligners["ajpolar"] = factory<ExternalAJPolarAlignPtr>();
  // Register external alignment a vector filed to the external align class factory
  external_aligners["field"] = factory<ExternalFieldAlignPtr>();
  // Register active jamming nematic external alignment to the external align class factory
  external_aligners["ajnematic"] = factory<ExternalAJNematicAlignPtr>();
  // Register cell shape external alignment to the external align class factory
  external_aligners["cell_shape"] = factory<ExternalShapeAlignPtr>();
  // Register alignment to polymer tangent to the external align class factory
  external_aligners["tangent"] = factory<ExternalTangentAlignPtr>();
  // Register kenotaxis alignment class
  external_aligners["kenotaxis"] = factory<ExternalKenotaxisAlignPtr>();
  // Register radial alignment class
  external_aligners["radial"] = factory<ExternalRadialAlignPtr>();
  // Register piv alignment class
  external_aligners["piv"] = factory<ExternalPIVAlignPtr>();
}
