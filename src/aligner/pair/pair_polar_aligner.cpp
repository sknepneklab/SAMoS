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
 * \file pair_polar_aligner.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Jan-2014
 * \brief Implementation of PairPolarAlign class
 */ 

#include "pair_polar_aligner.hpp"

void PairPolarAlign::compute()
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double J = m_J;
  double rcut = m_rcut;
  
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_align_energy("polar",0.0);
    }
  }
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      if (m_has_pair_params)
      {
        int pi_t = pi.get_type() - 1, pj_t = pj.get_type() - 1;
        J = m_pair_params[pi_t][pj_t].J;
        rcut = m_pair_params[pi_t][pj_t].rcut;
      }
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
      if (r_sq <= rcut*rcut)
      {
        double tau_x = pi.ny*pj.nz - pi.nz*pj.ny;
        double tau_y = pi.nz*pj.nx - pi.nx*pj.nz;
        double tau_z = pi.nx*pj.ny - pi.ny*pj.nx;
        pi.tau_x +=  J*tau_x;
        pi.tau_y +=  J*tau_y;
        pi.tau_z +=  J*tau_z;
        pj.tau_x += -J*tau_x;
        pj.tau_y += -J*tau_y;
        pj.tau_z += -J*tau_z;
        double potential_energy = -2.0*J*(pi.nx*pj.nx + pi.ny*pj.ny + pi.nz*pj.nz);  // 2.0 needed since we only use half of the neighbour list
        m_potential_energy += potential_energy;
        if (m_system->compute_per_particle_energy())
        {
          pi.add_align_energy("polar",potential_energy);
          pj.add_align_energy("polar",potential_energy);
        }
      }
    }
  }
}