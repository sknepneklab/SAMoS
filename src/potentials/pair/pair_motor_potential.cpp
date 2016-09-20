/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015, 2016
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015, 2016
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file pair_motor_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-Sep-2016
 * \brief Implementation of PairMotorPotential class
 */ 

#include "pair_motor_potential.hpp"

//! \param dt time step sent by the integrator 
void PairMotorPotential::compute(double dt)
{
  int N = m_system->size();
  double alpha = m_alpha;
  double ai, aj;
  double force_factor;
  double phi_i = 1.0;  // phase in factor for particle i
  double phi_j = 1.0;  // phase in factor for particle j
  double phi = 1.0;    // phase in factor for pair interaction (see below)
  
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("motor",0.0);
    }
  }
  
  m_potential_energy = 0.0;
  double tot_pot = m_potential_energy;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_phase_in)
      phi_i = 0.5*(1.0 + m_val->get_val(static_cast<int>(pi.age/dt)));
    ai = pi.get_radius();
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      double n_dot_n = pi.nx*pj.nx + pi.ny*pj.ny + pi.nz*pj.nz;
      if (n_dot_n < 0.0)  // particles point in the "opposite" direction
      {
        if (m_phase_in)
        {
          phi_j = 0.5*(1.0 + m_val->get_val(static_cast<int>(pj.age/dt)));
          // Determine global phase in factor: particles start at 0.5 strength (both daugthers of a division replace the mother)
          // Except for the interaction between daugthers which starts at 0
          if (phi_i < 1.0 && phi_j < 1.0)
            phi = phi_i + phi_j - 1.0;
          else 
            phi = phi_i*phi_j;
        }
        alpha = m_pair_params[pi.get_type()-1][pj.get_type()-1].alpha;
        aj = pj.get_radius();
        double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
        m_system->apply_periodic(dx,dy,dz);
        double r_sq = dx*dx + dy*dy + dz*dz;
        double r = sqrt(r_sq);
        double ai_p_aj;
        if (!m_use_particle_radii)
          ai_p_aj = m_pair_params[pi.get_type()-1][pj.get_type()-1].a;
        else
          ai_p_aj = ai+aj;
        if (r < ai_p_aj)
        {
          force_factor = alpha*phi*n_dot_n;
          // Handle force
          pi.fx += force_factor*pi.nx;
          pi.fy += force_factor*pi.ny;
          pi.fz += force_factor*pi.nz;
          // This breaks 3d Newton's law!!!
          pj.fx += force_factor*pj.nx;
          pj.fy += force_factor*pj.ny;
          pj.fz += force_factor*pj.nz;
        }
      }
    }
  }
  m_potential_energy = tot_pot;
}
