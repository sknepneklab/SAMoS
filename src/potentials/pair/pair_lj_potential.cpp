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
 * \file pair_lj_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Implementation of PairLJPotential class
 */ 

#include "pair_lj_potential.hpp"

void PairLJPotential::compute(double dt)
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double sigma = m_sigma;
  double eps = m_eps;
  double rcut = m_rcut;
  double sigma_sq = sigma*sigma, rcut_sq = rcut*rcut;
  double alpha = 1.0;    // Phase in factor 
 
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("lj",0.0);
    }
  }

  // Reset total potential energy to zero
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_phase_in)
      alpha = m_val->get_val(static_cast<int>(pi.age/dt));
    double ai = pi.get_radius();
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      if (m_has_pair_params)
      {
        rcut = m_pair_params[pi.get_type()-1][pj.get_type()-1].rcut;
        rcut_sq = rcut*rcut;
      }
      double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
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
      if (r_sq <= rcut_sq)
      {
        if (m_has_pair_params)
        {
          int pi_t = pi.get_type() - 1, pj_t = pj.get_type() - 1;
          sigma = m_pair_params[pi_t][pj_t].sigma;
          eps = m_pair_params[pi_t][pj_t].eps;
          sigma_sq = sigma*sigma;
        }
        if (m_use_particle_radii)
        {
          sigma = ai+pj.get_radius();
          sigma_sq = sigma*sigma;
        }
        double inv_r_sq = sigma_sq/r_sq;
        double inv_r_6  = inv_r_sq*inv_r_sq*inv_r_sq;
        // Handle potential 
        double potential_energy = 4.0*eps*alpha*inv_r_6*(inv_r_6 - 1.0);
        if (m_shifted)
        {
          double inv_r_cut_sq = sigma_sq/rcut_sq;
          double inv_r_cut_6 = inv_r_cut_sq*inv_r_cut_sq*inv_r_cut_sq;
          potential_energy -= 4.0 * eps * alpha * inv_r_cut_6 * (inv_r_cut_6 - 1.0);
        }
        m_potential_energy += potential_energy;
        // Handle force
        double force_factor = 48.0*eps*alpha*inv_r_6*(inv_r_6 - 0.5)*inv_r_sq;
        pi.fx += force_factor*dx;
        pi.fy += force_factor*dy;
        pi.fz += force_factor*dz;
        // Use 3d Newton's law
        pj.fx -= force_factor*dx;
        pj.fy -= force_factor*dy;
        pj.fz -= force_factor*dz;
        if (m_system->compute_per_particle_energy())
        {
          pi.add_pot_energy("lj",potential_energy);
          pj.add_pot_energy("lj",potential_energy);
        }
      }
    }
  }
}