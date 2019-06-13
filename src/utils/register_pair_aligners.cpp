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
 * \file regster_pair_aligners.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all pair aligners with the class factory.
*/

#include "register.hpp"

void register_pair_aligners(PairAlignerMap& pair_aligners)
{
  // Register polar aligner with the pairwise aligner class factory
  pair_aligners["polar"] = factory<PairPolarAlignPtr>();
  // Register nematic aligner with the pairwise aligner class factory
  pair_aligners["nematic"] = factory<PairNematicAlignPtr>();
  // Register Vicsek aligner with the pairwise aligner class factory
  pair_aligners["vicsek"] = factory<PairVicsekAlignPtr>();
  // Register velocity aligner with the pairwise aligner class factory
  pair_aligners["velocity"] = factory<PairVelocityAlignPtr>();
}
