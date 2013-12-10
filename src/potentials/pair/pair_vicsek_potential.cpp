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
 * \file pair_viscek_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 10-Dec-2013
 * \brief Implementation of PairVicsekPotential class
 */ 

#include "pair_vicsek_potential.hpp"

void PairVicsekPotential::compute()
{
  int N = m_system->size();
  
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      pi.tau_x += pj.vx;
      pi.tau_y += pj.vy;
      pi.tau_z += pj.vz;
      // Since we use 3d Newton's and only have half the the neighbour list we need handle the other part as well
      pj.tau_x += pi.vx;
      pj.tau_y += pi.vy;
      pj.tau_z += pi.vz;
    }
  }
  // Now we need to actually normalize all taus
  for  (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    double v = sqrt(p.tau_x*p.tau_x + p.tau_y*p.tau_y + p.tau_z*p.tau_z);
    if (v == 0.0) v = 1.0;
    double f = m_v0/v;
    p.tau_x *= f;
    p.tau_y *= f;
    p.tau_z *= f;
  }
  
}