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

//! Get cell to which given particle belongs
//! \param p Reference to the particle object
int CellList::get_cell_idx(const Particle& p)
{
  BoxPtr box = m_system->get_box();
  int i = static_cast<int>((p.x-box->xlo)/m_wx);
  int j = static_cast<int>((p.y-box->ylo)/m_wy);
  int k = static_cast<int>((p.z-box->zlo)/m_wz); 
  return m_ny*m_nz*i + m_nz*j + k;
}

//! Populate cell list
void CellList::populate()
{
  BoxPtr box = m_system->get_box();
  double lx = box->Lx, ly = box->Ly, lz = box->Lz;
  bool periodic = m_system->get_periodic();
  for (int i = 0; i < m_size; i++)
    m_cells[i].wipe();
  for (int i = 0; i < m_system->size(); i++)
  {
    Particle& p = m_system->get_particle(i);
    if (periodic)
    {
      if (p.x < box->xlo) p.x += lx;
      else if (p.x > box->xhi) p.x -= lx;
      if (p.y < box->ylo) p.y += ly;
      else if (p.y > box->yhi) p.y -= ly;
      if (p.z < box->zlo) p.z += lz;
      else if (p.z > box->zhi) p.z -= lz;
    }
    this->add_particle(p);
  }
  //m_msg->msg(Messenger::INFO,"Populated cell list.");
}
