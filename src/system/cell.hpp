/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

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