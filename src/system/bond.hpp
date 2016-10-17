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
 * \file bond.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Nov-2014
 * \brief Declaration of Bond data structure.
 */ 

#ifndef __BOND_HPP__
#define __BOND_HPP__

const int NUM_BOND_ATTRIB = 4;  //!< Number of bond attributes (columns in the bond file)

/*! Bond structure
 *  This is a basic data structure that hold information about bonds
 *  in the system. This is used to model chain-like molecules, such as
 *  polymers.
 *  Bond holds two indices, i and j which are id's of particles
 *  that are bonded. It also includes bond type, which is a positive
 *  integer.
 */
struct Bond
{
  //! Construct a bond object
  //! \param id Unique id of the bond
  //! \param type Bond type
  //! \param i Id of the first particle in the pair
  //! \param j Id of the second particle in the pair
  Bond(int id, int type, int i, int j) : id(id), type(type), i(i), j(j)  { }
  int id;     //!< Unique id of a bond
  int type;   //!< Bond type
  int i;      //!< Id of the first particle in the pair
  int j;      //!< Id of the second particle in the pair
};


#endif
