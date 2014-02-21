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
 *   (c) 2013
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file neighbour_list.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Oct-2013
 * \brief Definition of the NeighbourList class
*/

#include "neighbour_list.hpp"

/*! Build neighbour list using N^2 algorithm 
*/
void NeighbourList::build()
{
 if (m_use_cell_list) 
    this->build_cell();
  else
    this->build_nsq();
}

/* Do actual building. */

//! Builds neighbour list using N^2 (all pairs algorithm).
void NeighbourList::build_nsq()
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double cut = m_cut+m_pad;
  double cut2 = cut*cut;
  double d2;
  
  m_old_state.clear();
  
  for (int i = 0; i < N; i++)
    m_list[i].clear();
  
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    for (int j = i + 1; j < N; j++)
    {
      Particle& pj = m_system->get_particle(j);
      double dx = pi.x - pj.x;
      double dy = pi.y - pj.y;
      double dz = pi.z - pj.z;
      if (periodic)
      {
        if (dx > box->xhi) dx -= box->Lx;
        else if (dx < box->xlo) dx += box->Lx;
        if (dy > box->yhi) dy -= box->Ly;
        else if (dy < box->ylo) dy += box->Ly;
        if (dz > box->zhi) dz -= box->Lz;
        else if (dz < box->zlo) dz += box->Lz;
      }
      d2 = dx*dx + dy*dy + dz*dz;
      if (d2 < cut2)
        m_list[i].push_back(j);
    }
    m_old_state.push_back(PartPos(pi.x,pi.y,pi.z));
  }
  //m_msg->msg(Messenger::INFO, "Rebuilt neighbour list (N^2 algorithm).");
}

//! Build neighbour list using cell list
void NeighbourList::build_cell()
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double cut = m_cut+m_pad;
  double cut2 = cut*cut;
  double d2;
  
  m_old_state.clear();
  
  for (int i = 0; i < N; i++)
    m_list[i].clear();
  
  m_cell_list->populate();
  
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    int cell_idx = m_cell_list->get_cell_idx(pi);
    vector<int>& neigh_cells = m_cell_list->get_cell(cell_idx).get_neighbours();  // per design includes this cell as well
    for (vector<int>::iterator it = neigh_cells.begin(); it != neigh_cells.end(); it++)
    {
      Cell& c = m_cell_list->get_cell(*it);
      vector<int>& p_idx_vec = c.get_particles(); 
      for (unsigned int j = 0; j < p_idx_vec.size(); j++)
      {
        Particle& pj = m_system->get_particle(p_idx_vec[j]);
        if (pj.get_id() > pi.get_id())
        {
          double dx = pi.x - pj.x;
          double dy = pi.y - pj.y;
          double dz = pi.z - pj.z;
          if (periodic)
          {
            if (dx > box->xhi) dx -= box->Lx;
            else if (dx < box->xlo) dx += box->Lx;
            if (dy > box->yhi) dy -= box->Ly;
            else if (dy < box->ylo) dy += box->Ly;
            if (dz > box->zhi) dz -= box->Lz;
            else if (dz < box->zlo) dz += box->Lz;
          }
          d2 = dx*dx + dy*dy + dz*dz;
          if (d2 < cut2)
            m_list[i].push_back(pj.get_id());
        }
      }
    }
    m_old_state.push_back(PartPos(pi.x,pi.y,pi.z));
  }
  //m_msg->msg(Messenger::INFO, "Rebuilt neighbour list (using cell list).");
}
