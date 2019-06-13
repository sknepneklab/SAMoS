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
 * \file regster_integrators.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all integrators with the class factory.
*/

#include "register.hpp"

void register_integrators(IntegratorMap& integrators)
{
  // Register Brownian dynamics integrator with the integrators class factory
  integrators["brownian"] = factory<IntegratorBrownianPtr>();
  // Register Vicsek dynamics integrator with the integrators class factory
  integrators["vicsek"] = factory<IntegratorVicsekPtr>();
  // Register NVE integrator with the integrators class factory
  integrators["nve"] = factory<IntegratorNVEPtr>();
  // Register nematic integrator with the integrators class factory
  integrators["nematic"] = factory<IntegratorNematicPtr>();
  // Register actomyo integrator with the integrators class factory
  integrators["actomyo"] = factory<IntegratorActomyoPtr>();
  // Register Brownian rod dynamics integrator with the integrators class factory
  integrators["brownian_rod"] = factory<IntegratorBrownianRodPtr>();
  // Register Brownian dynamics integrator for particle position with the integrators class factory
  integrators["brownian_pos"] = factory<IntegratorBrownianPosPtr>();
  // Register Brownian dynamics integrator for alignment with the integrators class factory
  integrators["brownian_align"] = factory<IntegratorBrownianAlignPtr>();
  // Register Langevin dynamics integrator for particle positions with the integrators class factory
  integrators["langevin"] = factory<IntegratorLangevinPtr>();
  // Register FIRE minimisaton integrator with the integrators class factory
  integrators["fire"] = factory<IntegratorFIREPtr>();
  // Register Sepulveda minimisaton integrator with the integrators class factory
  integrators["sepulveda"] = factory<IntegratorSepulvedaPtr>();
}
