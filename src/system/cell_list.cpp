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
 * \file cell_list.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Jan-2014
 * \brief Implementation of CellList class members.
 */ 

#include "cell_list.hpp"

//! Construct cell list
//! \param sys Pointer to the system object
//! \param msg Pointer to the messenger object
//! \param cutoff cell size (currently all cells are cubic)
CellList::CellList(SystemPtr sys, MessengerPtr msg, double cutoff) : m_system(sys), m_msg(msg)
{
  m_nx = static_cast<int>(sys->get_box()->Lx/cutoff);
  m_ny = static_cast<int>(sys->get_box()->Ly/cutoff);
  m_nz = static_cast<int>(sys->get_box()->Lz/cutoff);
  m_size = m_nx*m_ny*m_nz;
  m_wx = sys->get_box()->Lx/m_nx;
  m_wy = sys->get_box()->Ly/m_ny;
  m_wz = sys->get_box()->Lz/m_nz;
  for (int i = 0; i < m_nx; i++)
    for (int j = 0; j < m_ny; j++)
      for (int k = 0; k < m_nz; k++)
      {
        int idx = m_ny*m_nz*i + m_nz*j + k;
        m_cells.push_back(Cell(idx));
        for (int ix = -1; ix <= 1; ix++)
          for (int iy = -1; iy <= 1; iy++)
            for (int iz = -1; iz <= 1; iz++)
            {
              int iix = i + ix, iiy = j + iy, iiz = k + iz;
              if (iix < 0) iix = m_nx - 1;
              else if (iix == m_nx) iix = 0;
              if (iiy < 0) iiy = m_ny - 1;
              else if (iiy == m_ny) iiy = 0;
              if (iiz < 0) iiz = m_nz - 1;
              else if (iiz == m_nz) iiz = 0;
              m_cells[idx].add_neighbour(m_ny*m_nz*iix + m_nz*iiy + iiz);
            }
      } 
  m_msg->msg(Messenger::INFO,"Created cell list with "+lexical_cast<string>(m_size)+" cells.");
  m_msg->msg(Messenger::INFO,"Each cell has dimensions ("+lexical_cast<string>(m_wx)+","+lexical_cast<string>(m_wy)+","+lexical_cast<string>(m_wz)+").");
}

//! Populate cell list
void CellList::populate()
{
  for (int i = 0; i < m_size; i++)
    m_cells[i].wipe();
  for (int i = 0; i < m_system->size(); i++)
  {
    Particle& p = m_system->get_particle(i);
    this->add_particle(p);
  }
  m_msg->msg(Messenger::INFO,"Populated cell list.");
}