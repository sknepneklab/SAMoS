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
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file pair_coulomb_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Implementation of PairCoulombPotential class
 */ 

#include "pair_coulomb_potential.hpp"

void PairCoulombPotential::compute(double dt)
{
  int N = m_system->size();
  double sigma = m_sigma;
  double alpha = m_alpha;
  double sigma_sq = sigma*sigma;
  double phase_fact = 1.0;   // Phase in factor
  
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
    if (m_phase_in)
      phase_fact = m_val->get_val(static_cast<int>(pi.age/dt));
    for (int j = i+1; j < N; j++)
    {
      Particle& pj = m_system->get_particle(j);
      if (m_has_pair_params)
      {
        int pi_t = pi.get_type() - 1, pj_t = pj.get_type() - 1;
        alpha = m_pair_params[pi_t][pj_t].alpha;
        sigma = m_pair_params[pi_t][pj_t].sigma;
        sigma_sq = sigma*sigma;
      }
      double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
      m_system->apply_periodic(dx,dy,dz);
      double r_sq = dx*dx + dy*dy + dz*dz;
      double r = sqrt(r_sq);
      double inv_r_sq = sigma_sq/r_sq;
      double inv_r_6  = inv_r_sq*inv_r_sq*inv_r_sq;
      // Handle potential 
      double potential_energy = phase_fact*(alpha/r + 4.0*fabs(alpha)*inv_r_6*inv_r_6);
      m_potential_energy += potential_energy;
      // Handle force
      double r_3 = r*r_sq;
      double force_factor = phase_fact*(alpha/r_3 + 48.0*fabs(alpha)*inv_r_6*inv_r_6*inv_r_sq);
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