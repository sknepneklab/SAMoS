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
 *   (c) 2013, 2014, 2015
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file pair_rod_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Implementation of PairRodPotential class
 */ 

#include "pair_rod_potential.hpp"

//! \param dt time step sent by the integrator 
void PairRodPotential::compute(double dt)
{
  int N = m_system->size();
  double k = m_k;
  double ai, aj;
  double li, lj;
  double si, sj;
  double dx, dy, dz;
  double force_factor;
  double alpha = 1.0;  // phase in factor
  
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("rod",0.0);
    }
  }
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_phase_in)
      alpha = m_val->get_val(static_cast<int>(pi.age/dt));
    ai = pi.get_radius();
    li = pi.get_length();
    double ni_x = pi.nx, ni_y = pi.ny, ni_z = pi.nz; 
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      if (m_has_pair_params)
        k = m_pair_params[pi.get_type()-1][pj.get_type()-1].k;
      aj = pj.get_radius();
      lj = pj.get_length();
      double nj_x = pj.nx, nj_y = pj.ny, nj_z = pj.nz;
      double dx_cm = pj.x - pi.x, dy_cm = pj.y - pi.y, dz_cm = pj.z - pi.z;
      m_system->apply_periodic(dx_cm,dy_cm,dz_cm);
      
      double ni_dot_nj = ni_x*nj_x + ni_y*nj_y + ni_z*nj_z; 
      double drcm_dot_ni = dx_cm*ni_x + dy_cm*ni_y + dz_cm*ni_z;
      double drcm_dot_nj = dx_cm*nj_x + dy_cm*nj_y + dz_cm*nj_z;
      
      if (1.0 - fabs(ni_dot_nj) == 1e-4) // rods a parallel
      {
        si =  li/(li+lj)*drcm_dot_ni;
        sj = -lj/(li+lj)*drcm_dot_nj;
      }
      else    // rods are not parallel
      {
        double denom = 1.0/(1.0-ni_dot_nj*ni_dot_nj);
        si = -denom*(drcm_dot_nj*ni_dot_nj-drcm_dot_ni)/li;
        sj =  denom*(drcm_dot_ni*ni_dot_nj-drcm_dot_nj)/lj;
      }
      if (si > 0.5) si = 0.5;
      else if (si < -0.5) si = -0.5;
      if (sj > 0.5) sj = 0.5;
      else if (sj < -0.5) sj = -0.5;
      dx = dx_cm + sj*lj*nj_x - si*li*ni_x;
      dy = dy_cm + sj*lj*nj_y - si*li*ni_y;
      dz = dz_cm + sj*lj*nj_z - si*li*ni_z;
      m_system->apply_periodic(dx,dy,dz);
      
      double r_sq = dx*dx + dy*dy + dz*dz;
      double r = sqrt(r_sq);
      double ai_p_aj = ai+aj;
            
      if (r < ai_p_aj)
      {
        
        // Handle potential 
        double diff = ai_p_aj - r;
        double pot_eng = 0.5*k*alpha*diff*diff;
        m_potential_energy += pot_eng;
        // Handle force
        if (r > 0.0) force_factor = k*alpha*diff/r;
        else force_factor = k*alpha*diff;
        double fx = force_factor*dx, fy = force_factor*dy, fz = force_factor*dz;
        pi.fx -= fx;
        pi.fy -= fy;
        pi.fz -= fz;
        // Use 3d Newton's law
        pj.fx += fx;
        pj.fy += fy;
        pj.fz += fz;
        // handle torques
        pi.tau_x -= si*(ni_y*fz - ni_z*fy);
        pi.tau_y -= si*(ni_z*fx - ni_x*fz);
        pi.tau_z -= si*(ni_x*fy - ni_y*fx);
        pj.tau_x += sj*(nj_y*fz - nj_z*fy);
        pj.tau_y += sj*(nj_z*fx - nj_x*fz);
        pj.tau_z += sj*(nj_x*fy - nj_y*fx);
        if (m_system->compute_per_particle_energy())
        {
          pi.add_pot_energy("rod",pot_eng);
          pj.add_pot_energy("rod",pot_eng);
        }
      }
    }
  }
}