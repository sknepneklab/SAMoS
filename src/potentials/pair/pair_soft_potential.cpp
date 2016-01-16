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
 * \file pair_soft_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Implementation of PairSoftPotential class
 */ 

#include "pair_soft_potential.hpp"

//! \param dt time step sent by the integrator 
void PairSoftPotential::compute(double dt)
{
  int N = m_system->size();
  double k = m_k;
  double ai, aj;
  double force_factor;
  double alpha_i = 1.0;  // phase in factor for particle i
  double alpha_j = 1.0;  // phase in factor for particle j
  double alpha = 1.0;    // phase in factor for pair interaction (see below)
  
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("soft",0.0);
    }
  }
  
  m_potential_energy = 0.0;
  double tot_pot = m_potential_energy;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_phase_in)
      alpha_i = 0.5*(1.0 + m_val->get_val(static_cast<int>(pi.age/dt)));
    ai = pi.get_radius();
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      if (m_phase_in)
      {
        alpha_j = 0.5*(1.0 + m_val->get_val(static_cast<int>(pj.age/dt)));
        // Determine global phase in factor: particles start at 0.5 strength (both daugthers of a division replace the mother)
        // Except for the interaction between daugthers which starts at 0
        if (alpha_i < 1.0 && alpha_j < 1.0)
	        alpha = alpha_i + alpha_j - 1.0;
	      else 
	        alpha = alpha_i*alpha_j;
      }
      k = m_pair_params[pi.get_type()-1][pj.get_type()-1].k;
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
        // Handle potential 
        double diff = ai_p_aj - r;
        double pot_eng = 0.5*k*alpha*diff*diff;
        //m_potential_energy += pot_eng;
        tot_pot += pot_eng;
        // Handle force
        if (r > 0.0) force_factor = k*alpha*diff/r;
        else force_factor = k*diff;
        pi.fx -= force_factor*dx;
        pi.fy -= force_factor*dy;
        pi.fz -= force_factor*dz;
        // Use 3d Newton's law
        pj.fx += force_factor*dx;
        pj.fy += force_factor*dy;
        pj.fz += force_factor*dz;
        if (m_system->compute_per_particle_energy())
        {
          pi.add_pot_energy("soft",pot_eng);
          pj.add_pot_energy("soft",pot_eng);
        }
      }
    }
  }
  m_potential_energy = tot_pot;
}