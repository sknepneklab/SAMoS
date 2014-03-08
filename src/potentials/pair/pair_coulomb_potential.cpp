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
 * \file pair_coulomb_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Implementation of PairCoulombPotential class
 */ 

#include "pair_coulomb_potential.hpp"

void PairCoulombPotential::compute()
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double sigma = m_sigma;
  double alpha = m_alpha;
  double sigma_sq = sigma*sigma;
  
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("coulomb",0.0);
    }
  }
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    for (int j = i+1; j < N; j++)
    {
      Particle& pj = m_system->get_particle(j);
      if (m_has_pair_params)
      {
        alpha = m_pair_params[make_pair(pi.get_type(),pj.get_type())]["alpha"];
        sigma = m_pair_params[make_pair(pi.get_type(),pj.get_type())]["sigma"];
        sigma_sq = sigma*sigma;
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
      double r = sqrt(r_sq);
      double inv_r_sq = sigma_sq/r_sq;
      double inv_r_6  = inv_r_sq*inv_r_sq*inv_r_sq;
      // Handle potential 
      double potential_energy = alpha/r + 4.0*fabs(alpha)*inv_r_6*inv_r_6;
      m_potential_energy += potential_energy;
      // Handle force
      double r_3 = r*r_sq;
      double force_factor = alpha/r_3 + 48.0*fabs(alpha)*inv_r_6*inv_r_6*inv_r_sq;
      pi.fx -= force_factor*dx;
      pi.fy -= force_factor*dy;
      pi.fz -= force_factor*dz;
      // Use 3d Newton's law
      pj.fx += force_factor*dx;
      pj.fy += force_factor*dy;
      pj.fz += force_factor*dz;
      if (m_system->compute_per_particle_energy())
      {
        pi.add_pot_energy("coulomb",potential_energy);
        pj.add_pot_energy("coulomb",potential_energy);
      }
    }
  }
}