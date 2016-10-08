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
 * \file pair_boundary_bending_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Apr-2016
 * \brief Implementation of PairBoundaryBendingPotential class
 */ 

#include "pair_boundary_bending_potential.hpp"


//! \param dt time step sent by the integrator 
void PairBoundaryBendingPotential::compute(double dt)
{
  double kappa = m_kappa;
  double theta0 = m_theta0;
  int N = m_system->size();
  Mesh& mesh = m_system->get_mesh();
  
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("boundary_bending",0.0);
    }
  }
  
  m_potential_energy = 0.0;
  for (int v = 0; v < mesh.size(); v++)
  {
    Vertex& V = mesh.get_vertices()[v];
    if (V.boundary && V.attached)
    {
      Vertex& Vm = mesh.get_vertices()[V.neigh[0]];
      Vertex& Vp = mesh.get_vertices()[V.neigh[V.n_edges-1]];
      Particle& p  = m_system->get_particle(V.id);
      Particle& pm = m_system->get_particle(Vm.id);
      Particle& pp = m_system->get_particle(Vp.id);
      if (m_has_part_params)
      {
        kappa = m_particle_params[p.get_type() - 1].kappa;
        theta0 = m_particle_params[p.get_type() - 1].theta0;
      }
      Vector3d r_pm_p = Vm.r - V.r;
      Vector3d r_pp_p = Vp.r - V.r;
      double r_pm_p_len = r_pm_p.len();
      double r_pp_p_len = r_pp_p.len();
      double denom = 1.0/(r_pm_p_len*r_pp_p_len);
      double denom_2 = denom*denom;
      double vec_dot = dot(r_pm_p,r_pp_p);
      double ratio = denom*vec_dot;
      double theta = acos(ratio);
      double fact = 0;
      if (1.0 - ratio*ratio > 1e-7 )
        fact = kappa*(theta-theta0)/sqrt(1-ratio*ratio);
      
      Vector3d fm = fact*(denom*r_pp_p - denom_2*(vec_dot*r_pp_p_len)*r_pm_p.unit());
      Vector3d fp = fact*(denom*r_pm_p - denom_2*(vec_dot*r_pm_p_len)*r_pp_p.unit());
      Vector3d f = -(fm+fp);
      
      pm.fx += fm.x;
      pm.fy += fm.y;
      pm.fz += fm.z;
      
      pp.fx += fp.x;
      pp.fy += fp.y;
      pp.fz += fp.z;
      
      p.fx += f.x;
      p.fy += f.y;
      p.fz += f.z;
      
    }
  }
}
