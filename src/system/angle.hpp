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
 * \file angle.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Nov-2014
 * \brief Declaration of Angle data structure.
 */ 

#ifndef __ANGLE_HPP__
#define __ANGLE_HPP__


const int NUM_ANGLE_ATTRIB = 5;  //!< Number of angle attributes (columns in the angle file)

/*! Angle structure
 *  This is a basic data structure that hold information about angle
 *  in the system. This is used to model bending rigidity 
 *  in chain-like molecules, such as polymers.
 *  Angle holds three indices, i, j and k which are id's of particles
 *  that are bonded. It also includes bond type, which is a positive
 *  integer.
 */
struct Angle
{
  //! Construct a angle object
  //! \param id Unique id of the bond
  //! \param i Id of the first (left) particle 
  //! \param j Id of the second (middle) particle 
  //! \param k If of the third (right) particle
  //! \param type Bond type
  Angle(int id, int type, int i, int j, int k) : id(id), type(type), i(i), j(j), k(k) { }
  int id;     //!< Unique id of a bond
  int type;   //!< Angle type
  int i;      //!< Id of the first (left) particle
  int j;      //!< Id of the second (middle) particle
  int k;      //!< If of the third (right) particle 
};


#endif
