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
 * \file cell.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Jan-2014
 * \brief Declaration of Cell class.
 */ 

#ifndef __CELL_HPP__
#define __CELL_HPP__

#include <vector>

using std::vector;

/*! This class handles single element (cell) of the cell list 
*/
class Cell
{
public:
  
  //! Construct a single empty cell
  //! \param idx cell index
  Cell(int idx) : m_id(idx) { }
  
  //! Adds particle to the cell
  //! \param idx particle index
  void add_particle(int idx)
  {
    m_particles.push_back(idx);
  }
  
  //! Returns vector with all particles in the cell
  vector<int>& get_particles()
  {
    return m_particles;
  }
  
  //! Wipes particles vector (used for cell list rebuilds)
  void wipe()
  {
    m_particles.clear();
  }
  
  //! Adds index of the neighbouring cell
  //! \param idx particle index
  void add_neighbour(int idx)
  {
    m_neigh.push_back(idx);
  }
  
  //! Returns list of all neighbours
  vector<int>& get_neighbours() 
  {
    return m_neigh;
  }
  
private:
  
  int m_id;                     //!< Unique id (index) of the cell
  vector<int> m_particles;      //!< Stores indices of particles belonging to this cell
  vector<int> m_neigh;          //!< Stores indices of neighbouring cells
  
};


#endif
