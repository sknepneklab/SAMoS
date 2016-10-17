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
 * \file pair_vertex_particle_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Implementation of PairVertexParticlePotential class
 */ 

#include "pair_vertex_particle_potential.hpp"

static int prev(int i, int N)
{
  if (i == 0) return N-1;
  return i-1;
}

static int next(int i, int N)
{
  if (i == N-1) return 0;
  return i+1;
}

//! \param dt time step sent by the integrator 
void PairVertexParticlePotential::compute(double dt)
{
  int N = m_system->size();
  double K = m_K;
  double gamma = m_gamma;
  double lambda = m_lambda;
  double alpha = 1.0;  // phase in factor
  double pot_eng = 0.0;
  
  if (m_mesh_update_steps > 0)
    if (m_system->get_step() % m_mesh_update_steps == 0)
      m_nlist->build_mesh();
  
  Mesh& mesh = m_system->get_mesh();
  
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("vp",0.0);
    }
  }
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    Vertex& vi = mesh.get_vertices()[i];
    Vector3d Nvec = Vector3d(pi.Nx, pi.Ny, pi.Nz);
    if (m_phase_in)
      alpha = m_val->get_val(static_cast<int>(pi.age/dt));
    // First handle the vertex itself
    if (pi.in_tissue && (m_include_boundary || !vi.boundary))
    {
      if (m_has_part_params) 
      {
        K  = m_particle_params[vi.type-1].K;
        gamma = m_particle_params[vi.type-1].gamma;
        lambda = m_particle_params[vi.type-1].lambda;
      }
      double dA = 0.0; 
      dA = (vi.area - pi.A0);
      double area_term = 0.5*K*dA;
      double perim_term = gamma*vi.perim;
      pot_eng = 0.5*(K*dA*dA+gamma*vi.perim*vi.perim);
      Vector3d area_vec(0.0,0.0,0.0);
      Vector3d perim_vec(0.0,0.0,0.0);
      Vector3d con_vec(0.0,0.0,0.0);
      for (int f = 0; f < vi.n_faces; f++)
      {
        int fid = vi.dual[f];
        int fid_m = vi.dual[prev(f,vi.n_faces)];
        int fid_p = vi.dual[next(f,vi.n_faces)];
        
        Face& f_nu_m = mesh.get_faces()[fid_m];
        Face& f_nu   = mesh.get_faces()[fid];
        Face& f_nu_p = mesh.get_faces()[fid_p];
         
        Vector3d& r_nu_m = f_nu_m.rc; 
        Vector3d& r_nu   = f_nu.rc; 
        Vector3d& r_nu_p = f_nu_p.rc; 
        
        Vector3d cross_prod_1(0.0,0.0,0.0);
        if (!(f_nu.is_hole || f_nu_p.is_hole)) cross_prod_1 = (cross(r_nu_p, Nvec))*f_nu.get_jacobian(i);
        if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_1 = cross_prod_1 - (cross(r_nu_m, Nvec))*f_nu.get_jacobian(i);
        area_vec = area_vec + cross_prod_1;  
        
        Vector3d cross_prod_2(0.0,0.0,0.0);
        if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_2 = ((r_nu - r_nu_m).unit())*f_nu.get_jacobian(i);
        if (!(f_nu_p.is_hole || f_nu.is_hole)) cross_prod_2 -= (((r_nu_p - r_nu).unit())*f_nu.get_jacobian(i));
        perim_vec = perim_vec + cross_prod_2; 
        
        // Adding force contributions from the dual edges associated with face f and face next(f,vi.n_faces)
        if (m_has_pair_params)
        {
          Vertex& vn = mesh.get_vertices()[vi.dual_neighbour_map[f]];
          lambda = m_pair_params[vi.type-1][vn.type-1].lambda;
        } 
        Vector3d cross_prod_3(0.0,0.0,0.0);
        if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_3 = lambda*(((r_nu - r_nu_m).unit())*f_nu.get_jacobian(i));
        con_vec = con_vec + cross_prod_3; 
        
        // Adding force contributions from the dual edges associated with face f and face prev(f,vi.n_faces)
        if (m_has_pair_params)
        {
          Vertex& vn = mesh.get_vertices()[vi.dual_neighbour_map[next(f,vi.n_faces)]];
          lambda = m_pair_params[vi.type-1][vn.type-1].lambda;
        } 
        Vector3d cross_prod_4(0.0,0.0,0.0);
        if (!(f_nu_p.is_hole || f_nu.is_hole)) cross_prod_4 = lambda*(((r_nu_p - r_nu).unit())*f_nu.get_jacobian(i));
        con_vec = con_vec - cross_prod_4; 
        
        pot_eng += lambda*(r_nu - r_nu_m).len();
      }
      // area term
      double area_fact = -alpha*area_term;
      pi.fx += area_fact*area_vec.x;
      pi.fy += area_fact*area_vec.y;
      pi.fz += area_fact*area_vec.z;
      // perimeter term 
      double perim_fact = -alpha*perim_term;
      pi.fx += perim_fact*perim_vec.x;
      pi.fy += perim_fact*perim_vec.y;
      pi.fz += perim_fact*perim_vec.z;
      // Add contractile term
      pi.fx -= alpha*con_vec.x;
      pi.fy -= alpha*con_vec.y;
      pi.fz -= alpha*con_vec.z;
      // add potential energy
      m_potential_energy += pot_eng;
    }
    // Now check neighbours
    for (int j = 0; j < vi.n_edges; j++)
    {
      Particle& pj = m_system->get_particle(vi.neigh[j]);
      Vertex& vj = mesh.get_vertices()[vi.neigh[j]];
      if (pj.in_tissue && (m_include_boundary || !vj.boundary))  
      {
        if (m_has_part_params) 
        {
          K  = m_particle_params[vj.type-1].K;
          gamma = m_particle_params[vj.type-1].gamma;
          lambda = m_particle_params[vj.type-1].lambda;
        }
        double dA = 0.0;
        dA = (vj.area - pj.A0);
        double area_term = 0.5*K*dA;
        double perim_term = gamma*vj.perim;
        Vector3d area_vec(0.0,0.0,0.0);
        Vector3d perim_vec(0.0,0.0,0.0);
        Vector3d con_vec(0.0,0.0,0.0);
        Vector3d Nvec = Vector3d(pj.Nx, pj.Ny, pj.Nz);
        
        for (int f = 0; f < vj.n_faces; f++)
        {
          int fid = vj.dual[f];
          Face& f_nu   = mesh.get_faces()[fid];
          if (f_nu.has_vertex(i))
          {
            int fid_m = vj.dual[prev(f,vj.n_faces)];
            int fid_p = vj.dual[next(f,vj.n_faces)];
            Face& f_nu_m = mesh.get_faces()[fid_m];
            Face& f_nu_p = mesh.get_faces()[fid_p];
            Vector3d& r_nu_m = f_nu_m.rc; 
            Vector3d& r_nu   = f_nu.rc; 
            Vector3d& r_nu_p = f_nu_p.rc; 

            Vector3d cross_prod_1(0.0,0.0,0.0);
            if (!(f_nu.is_hole || f_nu_p.is_hole)) cross_prod_1 = cross(r_nu_p, Nvec)*f_nu.get_jacobian(i);
            if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_1 = cross_prod_1 - cross(r_nu_m, Nvec)*f_nu.get_jacobian(i);
            area_vec = area_vec + cross_prod_1; 
            
            Vector3d cross_prod_2(0.0,0.0,0.0);
            if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_2 = ((r_nu - r_nu_m).unit())*f_nu.get_jacobian(i);
            if (!(f_nu_p.is_hole || f_nu.is_hole)) cross_prod_2 -= ((r_nu_p - r_nu).unit())*f_nu.get_jacobian(i);
            perim_vec = perim_vec + cross_prod_2; 
            
            // Adding force contributions from the dual edges associated with face f and face next(f,vi.n_faces)
            if (m_has_pair_params)
            {
              Vertex& vn = mesh.get_vertices()[vj.dual_neighbour_map[f]];
              lambda = m_pair_params[vj.type-1][vn.type-1].lambda;
            } 
            Vector3d cross_prod_3(0.0,0.0,0.0);
            if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_3 = lambda*(((r_nu - r_nu_m).unit())*f_nu.get_jacobian(i));
            con_vec = con_vec + cross_prod_3; 

            // Adding force contributions from the dual edges associated with face f and face next(f,vi.n_faces)
            if (m_has_pair_params)
            {
              Vertex& vn = mesh.get_vertices()[vj.dual_neighbour_map[next(f,vj.n_faces)]];
              lambda = m_pair_params[vj.type-1][vn.type-1].lambda;
            } 
            Vector3d cross_prod_4(0.0,0.0,0.0);
            if (!(f_nu_p.is_hole || f_nu.is_hole)) cross_prod_4 = lambda*(((r_nu_p - r_nu).unit())*f_nu.get_jacobian(i));
            con_vec = con_vec - cross_prod_4; 
            
          }
        }
        // area term
        double area_fact = -alpha*area_term;
        pi.fx += area_fact*area_vec.x;
        pi.fy += area_fact*area_vec.y;
        pi.fz += area_fact*area_vec.z;
        // perimeter term
        double perim_fact = -alpha*perim_term;
        pi.fx += perim_fact*perim_vec.x;
        pi.fy += perim_fact*perim_vec.y;
        pi.fz += perim_fact*perim_vec.z;
        // Add contractile term
        pi.fx -= alpha*con_vec.x;
        pi.fy -= alpha*con_vec.y;
        pi.fz -= alpha*con_vec.z;
        if (m_compute_stress)
        {
          if (!vi.boundary && vi.area > 0)
          {
            double inv_area = 1.0/vi.area;
            pi.s_xx *= inv_area;  pi.s_xy *= inv_area; pi.s_xz *= inv_area;
            pi.s_yx *= inv_area;  pi.s_yy *= inv_area; pi.s_yz *= inv_area;
            pi.s_zx *= inv_area;  pi.s_zy *= inv_area; pi.s_zz *= inv_area;
          }
        }
      }
    }
    // potential energy
    if (m_system->compute_per_particle_energy())
      pi.add_pot_energy("vp",pot_eng);
  }
}
