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
  int N = m_system->size();
  double cut2 = m_cut*m_cut;
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
      d2 = dx*dx + dy*dy + dz*dz;
      if (d2 < cut2)
        m_list[i].push_back(j);
    }
    m_old_state.push_back(PartPos(pi.x,pi.y,pi.z));
  }
  m_msg->msg(Messenger::INFO, "Rebuilt neighbour list.");
}