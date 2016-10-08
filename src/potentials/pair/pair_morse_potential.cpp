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
 * \file pair_morse_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Oct-2014
 * \brief Implementation of PairMorsePotential class
 */ 

#include "pair_morse_potential.hpp"

void PairMorsePotential::compute(double dt)
{
  int N = m_system->size();
  double D = m_D;
  double a = m_a;
  double re = m_re;
  double rcut = m_rcut;
  double rcut_sq = rcut*rcut;
  double alpha_i = 1.0;  // phase in factor for particle i
  double alpha_j = 1.0;  // phase in factor for particle j
  double alpha = 1.0;    // phase in factor for pair interaction (see below)
 
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("morse",0.0);
    }
  }

  // Reset total potential energy to zero
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_phase_in)
      alpha_i = 0.5*(1.0 + m_val->get_val(static_cast<int>(pi.age/dt)));
    double ai = pi.get_radius();
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
      if (m_has_pair_params)
      {
        rcut = m_pair_params[pi.get_type()-1][pj.get_type()-1].rcut;
        rcut_sq = rcut*rcut;
      }
      double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
      m_system->apply_periodic(dx,dy,dz);
      double r_sq = dx*dx + dy*dy + dz*dz;
      if (r_sq <= rcut_sq)
      {
        if (m_has_pair_params)
        {
          int pi_t = pi.get_type() - 1, pj_t = pj.get_type() - 1;
          D = m_pair_params[pi_t][pj_t].D;
          a = m_pair_params[pi_t][pj_t].a;
          re = m_pair_params[pi_t][pj_t].re;
        }
        if (m_use_particle_radii)
          re = ai + pj.get_radius();
        double r = sqrt(r_sq);
        double exp_fact = exp(-a*(r-re));
        double pot_fact = exp_fact - 1.0;
        // Handle potential 
        double potential_energy = D*alpha*(pot_fact*pot_fact-1.0);
        if (m_shifted)
        {
          double shift_fact = exp(-a*(rcut-re)) - 1.0;
          potential_energy -= D*alpha*(shift_fact*shift_fact-1.0);
        }
        m_potential_energy += potential_energy;
        // Handle force
        double force_factor = 2.0*D*a*alpha*exp_fact*pot_fact/r;
        pi.fx += force_factor*dx;
        pi.fy += force_factor*dy;
        pi.fz += force_factor*dz;
        // Use 3d Newton's law
        pj.fx -= force_factor*dx;
        pj.fy -= force_factor*dy;
        pj.fz -= force_factor*dz;
        if (m_system->compute_per_particle_energy())
        {
          pi.add_pot_energy("morse",potential_energy);
          pj.add_pot_energy("morse",potential_energy);
        }
      }
    }
  }
}
