/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

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
  int get_cell_idx(const Particle& p)
  {
    int i = static_cast<int>((p.x-m_system->get_box()->xlo)/m_wx);
    int j = static_cast<int>((p.y-m_system->get_box()->ylo)/m_wy);
    int k = static_cast<int>((p.z-m_system->get_box()->zlo)/m_wz); 
    return m_ny*m_nz*i + m_nz*j + k;
  }
  
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