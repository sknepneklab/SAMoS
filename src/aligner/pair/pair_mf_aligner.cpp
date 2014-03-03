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

#include "pair_mf_aligner.hpp"

void PairMFAlign::compute()
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double J = m_J;
  double rcut = m_rcut;
  
  // Reset MF vector B
  for  (int i = 0; i < N; i++)
  {
    m_B[i].x = 0.0;  m_B[i].y = 0.0;  m_B[i].z = 0.0; 
  }
  
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    vector<int>& neigh = m_nlist->get_neighbours(i);
    double n_neigh = static_cast<double>(neigh.size());
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      if (m_has_pair_params)
      {
        J = m_pair_params[make_pair(pi.get_type(),pj.get_type())]["J"];
        rcut = m_pair_params[make_pair(pi.get_type(),pj.get_type())]["rcut"];
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
        m_B[i].x += J/n_neigh*pj.nx;  m_B[i].y += J/n_neigh*pj.ny;  m_B[i].z += J/n_neigh*pj.nz;
        m_B[j].x += J/n_neigh*pi.nx;  m_B[j].y += J/n_neigh*pi.ny;  m_B[j].z += J/n_neigh*pi.nz;
      }
    }
  }
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    pi.tau_x +=  m_B[i].z*pi.ny - m_B[i].y*pi.nz;
    pi.tau_y += -m_B[i].z*pi.nx + m_B[i].x*pi.nz;
    pi.tau_z +=  m_B[i].y*pi.nx - m_B[i].x*pi.ny;  
    m_potential_energy += sqrt(pi.tau_x*pi.tau_x + pi.tau_y*pi.tau_y + pi.tau_z*pi.tau_z);  // Hm... Is this one OK???
    //std::cout << i << "  " << pi.tau_x << "  " << pi.tau_y << "  " << pi.tau_z << std::endl;
  }
}