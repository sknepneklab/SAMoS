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
 * \file cell_list.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Jan-2014
 * \brief Declaration of CellList class.
 */ 

#ifndef __CELL_LIST_HPP__
#define __CELL_LIST_HPP__

#include <vector>
#include <stdexcept>

#include "messenger.hpp"
#include "system.hpp"
#include "cell.hpp"

using std::vector;

/*! Handles cell list in the system. For simplicity all cell have 
 *  same dimensions (lx,ly,lz) and are given in 3 dimensions.
*/
class CellList
{
public:
  
  //! Construct cell list
  CellList(SystemPtr, MessengerPtr, double); 
  
  //! Get cell to which given particle belongs
  int get_cell_idx(const Particle&);
  
  //! Return total number of cells
  int get_size() { return m_size; }
  
  //! Return reference to a particular cell
  Cell& get_cell(int idx)  {  return m_cells[idx];  }
  
  //! Add particle to the appropriate cell
  //! \param p Reference to the particle
  void add_particle(const Particle& p)
  {
    int idx = this->get_cell_idx(p);
    m_cells[idx].add_particle(p.get_id());
  }
  
  //! Populates cell list
  void populate();
  
private:
  
  SystemPtr m_system;              //!< Pointer to the System object
  MessengerPtr m_msg;              //!< Handles messages sent to output
  vector<Cell> m_cells;            //!< Vector holding all cells in the system
  int m_size;                      //!< Cell list size (number of cells)
  double m_nx, m_ny, m_nz;         //!< Number of cell is x, y, z direction
  double m_wx, m_wy, m_wz;         //!< Cell width in the x, y, and z direction
  
};

typedef shared_ptr<CellList> CellListPtr;

#endif
