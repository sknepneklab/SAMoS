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
 * \file pair_soft_attractive_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-May-2015
 * \brief Implementation of PairSoftAttractivePotential class
 */ 

#include "pair_soft_attractive_potential.hpp"

//! \param dt time step sent by the integrator 
void PairSoftAttractivePotential::compute(double dt)
{
  int N = m_system->size();
  double k = m_k;
  double ai, aj;
  double force_factor = 0.0;
  double alpha_i = 1.0;  // phase in factor for particle i
  double alpha_j = 1.0;  // phase in factor for particle j
  double alpha = 1.0; // phase in factor for pair interaction (see below)
  
  for  (int i = 0; i < N; i++)
    {
	  Particle& p = m_system->get_particle(i);
	  if (m_system->compute_per_particle_energy())
	  {
		  p.set_pot_energy("soft_attractive",0.0);
	   }
	  // currently just hacking (for consistency with Daniel's code): Coordination number is interaction number
	  // Reset to 0 here and then recalculation in interaction below
	  p.coordination=0;
   }
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    // A note on phasing in: m_val, the value object, has been pre-set with the number of phase in time steps
    // It will give a linear interpolation between 0 and 1 based on the current age of the particle in time steps
    // Second note: I could set the interpolation between 0.5 and 1, however then pair_vertex potential would be dubious
    // Instead this is done manually here
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
      double fact = m_pair_params[pi.get_type()-1][pj.get_type()-1].fact;
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
      double r_turn = fact*ai_p_aj;
      double r_f = 2.0*r_turn - ai_p_aj;
      double diff;
      double pot_eng = 0.0;
      if (r < r_f)
      {
		    pi.coordination+=1;
		    pj.coordination+=1;
        if (r <= r_turn)
        {
          diff = ai_p_aj - r;
          pot_eng = 0.5*k*alpha*diff*diff;
          if (r > 0.0) force_factor = k*alpha*diff/r;
          else force_factor = k*diff;
        }
        else if (r > r_turn && r < r_f)
        {
          diff = r_f - r;
          pot_eng = -0.5*k*alpha*diff*diff;
          force_factor = -k*alpha*diff/r;
        }
        m_potential_energy += pot_eng;
        // Handle force
        pi.fx -= force_factor*dx;
        pi.fy -= force_factor*dy;
        pi.fz -= force_factor*dz;
        // Use 3d Newton's law
        pj.fx += force_factor*dx;
        pj.fy += force_factor*dy;
        pj.fz += force_factor*dz;
        if (m_system->compute_per_particle_energy())
        {
          pi.add_pot_energy("soft_attractive",pot_eng);
          pj.add_pot_energy("soft_attractive",pot_eng);
        }
      }
    }
  }
}
