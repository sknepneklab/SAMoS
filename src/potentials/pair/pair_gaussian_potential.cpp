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
 * \file pair_gaussian_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 27-Sept-2014
 * \brief Implementation of PairGaussianPotential class
 */ 

#include "pair_gaussian_potential.hpp"

void PairGaussianPotential::compute(double dt)
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double A = m_A;
  double B = m_B;
  double alpha = m_alpha;
  double beta = m_beta;
  double rA = m_rA;
  double rB = m_rB;
  double rcut = m_rcut;
  double rcut_sq = rcut*rcut;
  double phase_fact = 1.0;  // Phase in factor
    
   
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("gaussian",0.0);
    }
  }

  // Reset total potential energy to zero
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_phase_in)
      phase_fact = m_val->get_val(static_cast<int>(pi.age/dt));
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      double ai = pi.get_radius();
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
          A = m_pair_params[pi_t][pj_t].A;
          B = m_pair_params[pi_t][pj_t].B;
          alpha = m_pair_params[pi_t][pj_t].alpha;
          beta = m_pair_params[pi_t][pj_t].beta;
          rA = m_pair_params[pi_t][pj_t].rA;
          rB = m_pair_params[pi_t][pj_t].rB;
        }
        if (m_use_particle_radii)
          rB = ai + pj.get_radius();
        double r = sqrt(r_sq);
        // Handle potential 
        double r_m_rA = r - rA, r_m_rB = r - rB;
        double potential_energy = A*phase_fact*exp(-alpha*r_m_rA*r_m_rA)+B*exp(-beta*r_m_rB*r_m_rB);
        if (m_shifted)
        {
          double rcut_m_rA = rcut - rA, rcut_m_rB = rcut - rB;
          potential_energy -= A*phase_fact*exp(-alpha*rcut_m_rA*rcut_m_rA)+B*exp(-beta*rcut_m_rB*rcut_m_rB);
        }
        m_potential_energy += potential_energy;
        // Handle force
        double force_factor = (2.0/r)*phase_fact*(A*alpha*r_m_rA*exp(-alpha*r_m_rA*r_m_rA)+B*beta*r_m_rB*exp(-beta*r_m_rB*r_m_rB));
        pi.fx += force_factor*dx;
        pi.fy += force_factor*dy;
        pi.fz += force_factor*dz;
        // Use 3d Newton's law
        pj.fx -= force_factor*dx;
        pj.fy -= force_factor*dy;
        pj.fz -= force_factor*dz;
        if (m_system->compute_per_particle_energy())
        {
          pi.add_pot_energy("gaussian",potential_energy);
          pj.add_pot_energy("gaussian",potential_energy);
        }
      }
    }
  }
}