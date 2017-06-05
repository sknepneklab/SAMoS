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
  double beta = m_beta;
  double ai, aj;
  double force_factor;
  double phi_i = 1.0;  // phase in factor for particle i
  double phi_j = 1.0;  // phase in factor for particle j
  double phi = 1.0;    // phase in factor for pair interaction (see below)
  double activity;     // actual activity on each bead
  
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("motor",0.0);
    }
  }

  // Compute centres of mass of all filaments
  for (int mol = 0; mol < m_system->number_of_molecules(); mol++)
  {
    double xcm = 0.0, ycm = 0.0, zcm = 0.0;
    m_system->molecule_cm(mol,xcm,ycm,zcm);
    m_fil_cm[mol].x = xcm; m_fil_cm[mol].y = ycm; m_fil_cm[mol].z = zcm;
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
      if (pi.molecule != pj.molecule)  // Only particles on different molecules can interact
      {
        double n_dot_n = pi.nx*pj.nx + pi.ny*pj.ny + pi.nz*pj.nz; 
        if (m_phase_in)
        {
          phi_j = 0.5*(1.0 + m_val->get_val(static_cast<int>(pj.age/dt)));
          // Determine global phase in factor: particles start at 0.5 strength (both daughters of a division replace the mother)
          // Except for the interaction between daughters which starts at 0
          if (phi_i < 1.0 && phi_j < 1.0)
            phi = phi_i + phi_j - 1.0;
          else 
            phi = phi_i*phi_j;
        }
        alpha = m_pair_params[pi.get_type()-1][pj.get_type()-1].alpha;
        beta = m_pair_params[pi.get_type()-1][pj.get_type()-1].beta;
        if (n_dot_n < 0.0)  // filaments point in the "opposite" direction
          activity = alpha;
        else                // filaments point in the same direction
          activity = beta;
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
        // new code 
        /* comment out from here... */
       //  if (r < ai_p_aj)
//         {
//           if (n_dot_n <= 0.0)  // opposite direction
//           {
//             force_factor = alpha*phi*fabs(n_dot_n);
//             // Handle force
//             pi.fx += force_factor*pi.nx;
//             pi.fy += force_factor*pi.ny;
//             pi.fz += force_factor*pi.nz;
//             // This breaks 3d Newton's law!!! 
//             pj.fx += force_factor*pj.nx;
//             pj.fy += force_factor*pj.ny;
//             pj.fz += force_factor*pj.nz;
//           }
//           else                 // same direction 
//           {
//             // Note: This is inefficient. We are computing the same distance over and over again
//             double Xij = m_fil_cm[pi.molecule].x - m_fil_cm[pj.molecule].x; 
//             double Yij = m_fil_cm[pi.molecule].y - m_fil_cm[pj.molecule].y;
//             double Zij = m_fil_cm[pi.molecule].z - m_fil_cm[pj.molecule].z;
//             m_system->apply_periodic(Xij,Yij,Zij);
//             double ni_dot_Rij = pi.nx*Xij + pi.ny*Yij + pi.nz*Zij;
//             double nj_dot_Rij = pj.nx*Xij + pj.ny*Yij + pj.nz*Zij;
//             double sgn_i = (ni_dot_Rij > 0) ? 1.0 : -1.0;
//             double sgn_j = (nj_dot_Rij > 0) ? 1.0 : -1.0;
//             force_factor = beta*phi;
//             // Handle force
//             pi.fx += -sgn_i*force_factor*pi.nx;
//             pi.fy += -sgn_i*force_factor*pi.ny;
//             pi.fz += -sgn_i*force_factor*pi.nz;
//             // This breaks 3d Newton's law!!! 
//             pj.fx += sgn_j*force_factor*pj.nx;
//             pj.fy += sgn_j*force_factor*pj.ny;
//             pj.fz += sgn_j*force_factor*pj.nz;
//           }
//         }
        /* ... to here  */
        // old code 
        //uncomment from here...
        if (r < ai_p_aj)
        {
          force_factor = activity*phi*fabs(n_dot_n);
          // Handle force
          pi.fx += force_factor*pi.nx;
          pi.fy += force_factor*pi.ny;
          pi.fz += force_factor*pi.nz;
          // This breaks 3d Newton's law!!!
          // Variant with every two filaments pushing each other back (to force nematic behaviour): 
          // Make the force along - \hat{n} for one. How to choose consistently which one?
          // Trial run: since j < i always (see neighbour list), and indices are frozen in time, 
          // it should be enough to simply put a minus in front of the j one.
          force_factor = ((n_dot_n > 0.0) && m_allpairpush) ? -force_factor : force_factor; 
          pj.fx += force_factor*pj.nx;
          pj.fy += force_factor*pj.ny;
          pj.fz += force_factor*pj.nz;
        } 
        //to here */ 
      }
    }
  }
  m_potential_energy = tot_pot;
}
