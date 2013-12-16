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
 * \file pair_soft_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Implementation of PairSoftPotential class
 */ 

#include "pair_soft_potential.hpp"

void PairSoftPotential::compute()
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double k = m_k;
  double ai, aj;
  double force_factor;
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    ai = pi.get_radius();
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      if (m_has_pair_params)
        k = m_pair_params[make_pair(pi.get_type(),pj.get_type())]["k"];
      aj = pj.get_radius();
      double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
      if (periodic)
      {
        if (dx > box->xhi) dx -= box->Lx;
        else if (dx < box->xlo) dx += box->Lx;
        if (dy > box->yhi) dy -= box->Ly;
        else if (dy < box->ylo) dy += box->Ly;
        if (dz > box->zhi) dz -= box->Lz;
        else if (dz < box->zlo) dz += box->Lz;
      }
      double r_sq = dx*dx + dy*dy + dz*dz;
      double r = sqrt(r_sq);
      double ai_p_aj = ai + aj;
      if (r < ai_p_aj)
      {
        // Handle potential 
        double diff = ai_p_aj - r;
        m_potential_energy += 0.5*k*diff*diff;
        // Handle force
        if (r > 0.0) force_factor = k*diff/r;
        else force_factor = k*diff;
        pi.fx -= force_factor*dx;
        pi.fy -= force_factor*dy;
        pi.fz -= force_factor*dz;
        // Use 3d Newton's law
        pj.fx += force_factor*dx;
        pj.fy += force_factor*dy;
        pj.fz += force_factor*dz;
      }
    }
  }
}