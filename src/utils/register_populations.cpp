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
 * \file regster_populations.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all populations with the class factory.
*/

#include "register.hpp"

void register_populations(PopulationMap& populations)
{
  // Register random population control with the class factory
  populations["random"] = factory<PopulationRandomPtr>();
  // Register density population control with the class factory
  populations["density"] = factory<PopulationDensityPtr>();
  // Register growth population control with the class factory
  populations["grow"] = factory<PopulationGrowthPtr>();
  // Register elongation population control with the class factory
  populations["elongate"] = factory<PopulationElongationPtr>();
  // Register cell population control with the class factory
  populations["cell"] = factory<PopulationCellPtr>();
  // Register actomyosin population control with the class factory
  populations["actomyosin"] = factory<PopulationActomyosinPtr>();
  // Register actomyosin poisson population control with the class factory
  populations["actomyosin_poisson"] = factory<PopulationActomyosinPoissonPtr>();
  // Register actomyosin molecule population control with the class factory
  populations["actomyosin_molecule"] = factory<PopulationActomyosinMoleculePtr>();
  // Register actomyosin head population control with the class factory
  populations["actomyosin_head"] = factory<PopulationActomyosinHeadPtr>();
  // Register population regions control with the class factory
  populations["region"] = factory<PopulationRegionPtr>();
}
