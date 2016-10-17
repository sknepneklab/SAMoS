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
 * \file aligner.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Jan-2014
 * \brief Definition of Aligner class members.
*/ 

#include "aligner.hpp"

/*! Iterate over all pair and external aligners and compute 
 *  torques
 */
void Aligner::compute()
{
  //m_system->reset_torques();
  PairAlignType::iterator it_pair;
  ExternAlignType::iterator it_ext;
  
  for(it_pair = m_pair_align.begin(); it_pair != m_pair_align.end(); it_pair++)
    (*it_pair).second->compute();
  for(it_ext = m_external_align.begin(); it_ext != m_external_align.end(); it_ext++)
    (*it_ext).second->compute();
}
