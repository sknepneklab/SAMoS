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
 * \file pair_nematic_aligner.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 05-Mar-2014
 * \brief Implementation of PairNematicAlign class
 */ 

#include "pair_nematic_aligner.hpp"

void PairNematicAlign::compute()
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double J = 2.0*m_J;  // factor of 2 comes form the expansion of sin(2x) = 2sin(x)cos(x)
  double rcut = m_rcut;
  
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
        J = 2.0*m_pair_params[make_pair(pi.get_type(),pj.get_type())]["J"];
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
        double ni_dot_nj = pi.nx*pj.nx + pi.ny*pj.ny + pi.nz*pj.nz;
        double tau_x = pi.ny*pj.nz - pi.nz*pj.ny;
        double tau_y = pi.nz*pj.nx - pi.nx*pj.nz;
        double tau_z = pi.nx*pj.ny - pi.ny*pj.nx;
        pi.tau_x +=  J*ni_dot_nj*tau_x;
        pi.tau_y +=  J*ni_dot_nj*tau_y;
        pi.tau_z +=  J*ni_dot_nj*tau_z;
        pj.tau_x += -J*ni_dot_nj*tau_x;
        pj.tau_y += -J*ni_dot_nj*tau_y;
        pj.tau_z += -J*ni_dot_nj*tau_z;
        m_potential_energy += -2.0*J*(2.0*ni_dot_nj*ni_dot_nj - 1.0); // (cos(2x) = 2cos^2(x) - 1; factor 2.0 needed since we only use half of the neighbour list
      }
    }
  }
}