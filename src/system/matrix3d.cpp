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
 * \file vector3d.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 23-Mar-2016
 * \brief Definitions of Matrix3d class auxiliary functions
 */ 

#include "matrix3d.hpp"

//! Output vector components
std::ostream& operator<<(std::ostream& os, const Matrix3d& m)
{
  os << "[";
  os << "[" << m.M[0][0] << "," << m.M[0][1] << "," << m.M[0][2] << "],";
  os << "[" << m.M[1][0] << "," << m.M[1][1] << "," << m.M[1][2] << "],";
  os << "[" << m.M[2][0] << "," << m.M[2][1] << "," << m.M[2][2] << "]";
  os << "]";
  return os;
}
