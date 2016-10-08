/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

/*!
 * \file pair_rod_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 23-Mar-2015
 * \brief Implementation of PairRodPotential class
 */ 

#include "pair_rod_potential.hpp"

//! \param dt time step sent by the integrator 
void PairRodPotential::compute(double dt)
{
  int N = m_system->size();
  double k = m_k;
  double push = m_push;
  double ai, aj;
  double li2, lj2;
  double lambda, mu;
  double dx, dy, dz;
  double fx, fy, fz;
  double force_factor;
  double pot_eng;
  double alpha_i = 1.0;  // phase in factor for particle i
  double alpha_j = 1.0;  // phase in factor for particle j
  double alpha = 1.0;    // phase in factor for pair interaction (see below)
      
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
      alpha_i = 0.5*(1.0 + m_val->get_val(static_cast<int>(pi.age/dt)));
    ai = pi.get_radius();
    li2 = 0.5*pi.get_length();
    double ni_x = pi.nx, ni_y = pi.ny, ni_z = pi.nz; 
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      if (m_phase_in)
      {
        alpha_j = 0.5*(1.0 + m_val->get_val(static_cast<int>(pj.age/dt)));
        // Determine global phase in factor: particles start at 0.5 strength (both daugthers of a division replace the mother)
        // Except for the interaction between daugthers which starts at 0
        if (alpha_i<1.0 && alpha_j < 1.0)
	        alpha = alpha_i + alpha_j - 1.0;
	      else 
	        alpha = alpha_i*alpha_j;
      }
      if (m_has_pair_params)
      {
        k = m_pair_params[pi.get_type()-1][pj.get_type()-1].k;
        push = m_pair_params[pi.get_type()-1][pj.get_type()-1].push;
      }
      aj = pj.get_radius();
      lj2 = 0.5*pj.get_length();
      double nj_x = pj.nx, nj_y = pj.ny, nj_z = pj.nz;
      double dx_cm = pj.x - pi.x, dy_cm = pj.y - pi.y, dz_cm = pj.z - pi.z;
      m_system->apply_periodic(dx_cm,dy_cm,dz_cm);
      
      double ni_dot_nj = ni_x*nj_x + ni_y*nj_y + ni_z*nj_z; 
      double drcm_dot_ni = dx_cm*ni_x + dy_cm*ni_y + dz_cm*ni_z;
      double drcm_dot_nj = dx_cm*nj_x + dy_cm*nj_y + dz_cm*nj_z;
      double cc = 1.0 - ni_dot_nj*ni_dot_nj;
      
      
      if (cc <= 1e-10) // rods are nearly parallel
      {
        if (fabs(drcm_dot_ni) > 1e-10)
        {
          lambda = copysign(li2,drcm_dot_ni);
          mu = lambda*ni_dot_nj - drcm_dot_nj;
          if (fabs(mu) > lj2) mu = copysign(lj2,mu);
        }
        else
        {
          lambda = 0.0;
          mu = 0.0;
        }
//         lambda = 0.5*drcm_dot_ni;
//         mu = -0.5*drcm_dot_nj;
//         if (fabs(lambda) > li2) lambda = copysign(li2,lambda); 
//         if (fabs(mu) > lj2) mu = copysign(lj2,mu);
      }
      else
      {
        lambda = (drcm_dot_ni-ni_dot_nj*drcm_dot_nj)/cc;
        mu = (-drcm_dot_nj+ni_dot_nj*drcm_dot_ni)/cc;
        if (fabs(lambda) > li2 || fabs(mu) > lj2)
        {
          if (fabs(lambda) - li2 > fabs(mu) - lj2)
          {
            lambda = copysign(li2,lambda);
            mu = lambda*ni_dot_nj - drcm_dot_nj;
            if (fabs(mu) > lj2) mu = copysign(lj2,mu);
          }
          else
          {
            mu =  copysign(lj2,mu);
            lambda = mu*ni_dot_nj + drcm_dot_ni;
            if (fabs(lambda) > li2) lambda = copysign(li2,lambda);      
          }
        }
      }
           
      dx = dx_cm + mu*nj_x - lambda*ni_x;
      dy = dy_cm + mu*nj_y - lambda*ni_y;
      dz = dz_cm + mu*nj_z - lambda*ni_z;
      
      double r_sq = dx*dx + dy*dy + dz*dz;
      double r = sqrt(r_sq);
      double ai_p_aj = ai+aj;
     
      if (r < ai_p_aj)
      { 
        // Handle potential 
        double diff = ai_p_aj - r;
        pot_eng = 0.5*k*alpha*diff*diff;
        if (m_model == "hertz")
          pot_eng *= 0.8*sqrt(diff);
        m_potential_energy += pot_eng;
        // Handle force
        
        if (r >= SMALL_NUMBER)
        {
          force_factor = k*alpha*diff/r;
          if (m_model == "hertz")
            force_factor *= sqrt(diff);
          fx = force_factor*dx; fy = force_factor*dy; fz = force_factor*dz;
        }
        else
        {
          r_sq = dx_cm*dx_cm + dy_cm*dy_cm + dz_cm*dz_cm;
          r = sqrt(r_sq);
          if (r > 0.0)
            force_factor = k*alpha*push*diff/r;
          else
            force_factor = 0.0;
          if (m_model == "hertz")
            force_factor *= sqrt(diff);
          fx = force_factor*dx_cm; fy = force_factor*dy_cm ; fz = force_factor*dz_cm;
        }
        pi.fx -= fx;
        pi.fy -= fy;
        pi.fz -= fz;
        // Use 3d Newton's law
        pj.fx += fx;
        pj.fy += fy;
        pj.fz += fz;
        // handle torques
        pi.tau_x -= lambda*(ni_y*fz - ni_z*fy);
        pi.tau_y -= lambda*(ni_z*fx - ni_x*fz);
        pi.tau_z -= lambda*(ni_x*fy - ni_y*fx);
        pj.tau_x += mu*(nj_y*fz - nj_z*fy);
        pj.tau_y += mu*(nj_z*fx - nj_x*fz);
        pj.tau_z += mu*(nj_x*fy - nj_y*fx);
        if (m_system->compute_per_particle_energy())
        {
          pi.add_pot_energy("rod",pot_eng);
          pj.add_pot_energy("rod",pot_eng);
        }
      }
    }
  }
}
