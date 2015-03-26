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
 * \file pair_ljrod_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Mar-2015
 * \brief Implementation of PairLJRodPotential class
 */ 

#include "pair_ljrod_potential.hpp"

//! \param dt time step sent by the integrator 
void PairLJRodPotential::compute(double dt)
{
  int N = m_system->size();
  double sigma = m_sigma;
  double eps = m_eps;
  double rcut = m_rcut;
  double sigma_sq = sigma*sigma, rcut_sq = rcut*rcut;
  double li, lj;
  double si, sj;
  double dx, dy, dz;
  double fx, fy, fz;
  double force_factor;
  double potential_energy;
  double alpha = 1.0;  // phase in factor
  double k;
  double inv_core_sq, inv_core_6, lj_core_sq;
    
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("ljrod",0.0);
    }
  }
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_phase_in)
      alpha = m_val->get_val(static_cast<int>(pi.age/dt));
    li = pi.get_length();
    double ni_x = pi.nx, ni_y = pi.ny, ni_z = pi.nz; 
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      if (m_has_pair_params)
      {
        rcut = m_pair_params[pi.get_type()-1][pj.get_type()-1].rcut;
        rcut_sq = rcut*rcut;
      }
      lj = pj.get_length();
      double nj_x = pj.nx, nj_y = pj.ny, nj_z = pj.nz;
      double dx_cm = pj.x - pi.x, dy_cm = pj.y - pi.y, dz_cm = pj.z - pi.z;
      m_system->apply_periodic(dx_cm,dy_cm,dz_cm);
      
      double ni_dot_nj = ni_x*nj_x + ni_y*nj_y + ni_z*nj_z; 
      double drcm_dot_ni = dx_cm*ni_x + dy_cm*ni_y + dz_cm*ni_z;
      double drcm_dot_nj = dx_cm*nj_x + dy_cm*nj_y + dz_cm*nj_z;
      
      if (1.0 - fabs(ni_dot_nj) <= 1e-5) // rods are nearly parallel
      {
        si =  drcm_dot_ni/(li+lj);
        sj = -drcm_dot_nj/(li+lj);
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
      
      if (r_sq <= rcut_sq)
      { 
        if (m_has_pair_params)
        {
          int pi_t = pi.get_type() - 1, pj_t = pj.get_type() - 1;
          sigma = m_pair_params[pi_t][pj_t].sigma;
          eps = m_pair_params[pi_t][pj_t].eps;
          sigma_sq = sigma*sigma;
        }
        if (r_sq <= SMALL_NUMBER)
        {
          r_sq = dx_cm*dx_cm + dy_cm*dy_cm + dz_cm*dz_cm;
          dx = dx_cm;  dy = dy_cm;  dz = dz_cm;
          potential_energy = 4.0*eps;
          lj_core_sq = LJ_HARD_CORE_DISTANCE*LJ_HARD_CORE_DISTANCE;
          inv_core_sq = sigma_sq/lj_core_sq;
          inv_core_6 = inv_core_sq*inv_core_sq*inv_core_sq;
          k = 6.0*eps/(lj_core_sq*lj_core_sq)*inv_core_6*(2.0*inv_core_6 - 1.0);
          force_factor = 4.0*k*r_sq*sqrt(r_sq);   
          //cout << "Full overlap. k = " << k << endl;
        }
        else if (r_sq <= LJ_HARD_CORE_DISTANCE*LJ_HARD_CORE_DISTANCE*sigma_sq)
        {
          potential_energy = 4.0*eps;
          lj_core_sq = LJ_HARD_CORE_DISTANCE*LJ_HARD_CORE_DISTANCE;
          inv_core_sq = sigma_sq/lj_core_sq;
          inv_core_6 = inv_core_sq*inv_core_sq*inv_core_sq;
          k = 6.0*eps/(lj_core_sq*lj_core_sq)*inv_core_6*(2.0*inv_core_6 - 1.0);
          force_factor = 4.0*k*r_sq*sqrt(r_sq);
          //cout << "Too close. F = " << force_factor << endl;
        }
        else
        {
          double inv_r_sq = sigma_sq/r_sq;
          double inv_r_6  = inv_r_sq*inv_r_sq*inv_r_sq;
          // Handle potential 
          potential_energy = 4.0*eps*alpha*inv_r_6*(inv_r_6 - 1.0);
          if (m_shifted)
          {
            double inv_r_cut_sq = sigma_sq/rcut_sq;
            double inv_r_cut_6 = inv_r_cut_sq*inv_r_cut_sq*inv_r_cut_sq;
            potential_energy -= 4.0 * eps * alpha * inv_r_cut_6 * (inv_r_cut_6 - 1.0);
          }
          force_factor = 48.0*eps*alpha*inv_r_6*(inv_r_6 - 0.5)*inv_r_sq;
          //cout << "LJ part with F = " << force_factor << endl;
        }
        m_potential_energy += potential_energy;
        // Handle force
        fx = force_factor*dx;  fy = force_factor*dy;  fz = force_factor*dz;
        pi.fx -= fx;
        pi.fy -= fy;
        pi.fz -= fz;
        // Use 3d Newton's law
        pj.fx += fx;
        pj.fy += fy;
        pj.fz += fz;
        // handle torques
        pi.tau_x -= li*si*(ni_y*fz - ni_z*fy);
        pi.tau_y -= li*si*(ni_z*fx - ni_x*fz);
        pi.tau_z -= li*si*(ni_x*fy - ni_y*fx);
        pj.tau_x += lj*sj*(nj_y*fz - nj_z*fy);
        pj.tau_y += lj*sj*(nj_z*fx - nj_x*fz);
        pj.tau_z += lj*sj*(nj_x*fy - nj_y*fx);
        if (m_system->compute_per_particle_energy())
        {
          pi.add_pot_energy("ljrod",potential_energy);
          pj.add_pot_energy("ljrod",potential_energy);
        }
      }
    }
  }
}