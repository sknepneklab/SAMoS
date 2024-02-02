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
 * \file pair_abp_actreact_potential.hpp
 * \author Silke Henkes, silkehenkes@gmail.com
 * \date 01-February-2024
 * \brief Declaration of PairABPActReactPotential class 
 */ 

#include "pair_abp_actreact_potential.hpp"

//! \param dt time step sent by the integrator 
void PairABPActReactPotential::compute(double dt)
{
  int N = m_system->size();
  double poli, polj; // polarisation forces 
  double rinti, rintj; // interaction ranges multipliers
  double ai, aj;
  double alpha_i = 1.0;  // phase in factor for particle i
  double alpha_j = 1.0;  // phase in factor for particle j
  double alpha = 1.0; // phase in factor for pair interaction (see below)
  
  for  (int i = 0; i < N; i++)
    {
	  Particle& p = m_system->get_particle(i);
      // Not in fact computing energy for this non-conservative "potential"
	  if (m_system->compute_per_particle_energy())
	  {
		  p.set_pot_energy("abp_actreact",0.0);
	   }
   }
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    // type parameters
    std::cout << m_has_part_params << "type i " << pi.get_type() << endl;
    if (m_has_part_params)  
    {
      poli = m_particle_params[pi.get_type()-1].p;
      rinti = m_particle_params[pi.get_type()-1].r_int;
    }
    else 
    {
        poli = m_p;
        rinti = m_r_int;
    }
    
    // A note on phasing in: m_val, the value object, has been pre-set with the number of phase in time steps
    // It will give a linear interpolation between 0 and 1 based on the current age of the particle in time steps
    // Second note: I could set the interpolation between 0.5 and 1, however then pair_vertex potential would be dubious
    // Instead this is done manually here
    // The pair ABP inherits the whole phase-in mechanism from the other potentials, so that everything is as compatible as possible
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
      // data for particle 2
      aj = pj.get_radius();
      // type parameters
      if (m_has_part_params) 
      {
        polj = m_particle_params[pj.get_type()-1].p;
        rintj = m_particle_params[pj.get_type()-1].r_int;
      }
      else 
      {
        polj = m_p;
        rintj = m_r_int;
      }
      
      double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
      m_system->apply_periodic(dx,dy,dz);
      double r_sq = dx*dx + dy*dy + dz*dz;
      double r = sqrt(r_sq);

      double ai_p_aj;
      if (!m_use_particle_radii)
        ai_p_aj = 2.0;
      else
        ai_p_aj = ai+aj;

      double rintij = 0.5*(rinti+rintj)*ai_p_aj;
      double b;
      double fax, fay, faz;

      if (r < rintij)
      {
        b = (rintij-r)/rintij;  // force prefactor (positive and between 0 and 1)
        std::cout << "type j " << pj.get_type() << endl;
        std::cout << " pol i" << poli << ", rint i " << rinti << std::endl;
        std::cout << " pol j" << polj << ", rint j " << rintj << std::endl;
        std::cout << " beta " << b << "alpha " << alpha << std::endl;
        // action-reaction active force magnitude
        fax = 0.5*(poli*pi.nx-polj*pj.nx);
        fay = 0.5*(poli*pi.ny-polj*pj.ny);
        faz = 0.5*(poli*pi.nz-polj*pj.nz);
        // Handle force
        pi.fx += b*alpha*fax;
        pi.fy += b*alpha*fay;
        pi.fz += b*alpha*faz;
        // Use 3d Newton's law: this is an action-reaction and reciprocal active force
        pj.fx -= b*alpha*fax;
        pj.fy -= b*alpha*fay;
        pj.fz -= b*alpha*faz;
      }
      
    }
  }
}
