/* *************************************************************
 *  
 *   VertexParticle Active Mater on Surfaces (SAMoS)
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
    if (vi.attached && (m_include_boundary || !vi.boundary))
    {
      if (m_has_part_params) 
      {
        K  = m_particle_params[vi.type-1].K;
        gamma = m_particle_params[vi.type-1].gamma;
      }
      double dA = 0.0; 
      if (!vi.boundary) dA = (vi.area - pi.A0);
      else dA = (vi.area - mesh.angle_factor(i)*pi.A0);
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
        Vector3d r_nu_i = r_nu - vi.r;
        
        Vector3d cross_prod_1(0.0,0.0,0.0);
        if (!(f_nu.is_hole || f_nu_p.is_hole)) cross_prod_1 = (cross(r_nu_p, Nvec))*f_nu.get_jacobian(i);
        if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_1 = cross_prod_1 - (cross(r_nu_m, Nvec))*f_nu.get_jacobian(i);
        area_vec = area_vec + cross_prod_1;  
        
        if (vi.boundary)
        {
          if (f == 0)
            area_vec = area_vec + cross(r_nu,Nvec) - (cross(vi.r,Nvec))*f_nu.get_jacobian(i);
          else if (f == vi.n_faces - 2)
            area_vec = area_vec - cross(r_nu,Nvec) + (cross(vi.r,Nvec))*f_nu.get_jacobian(i);
        }
        Vector3d cross_prod_2(0.0,0.0,0.0);
        if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_2 = ((r_nu - r_nu_m).unit())*f_nu.get_jacobian(i);
        if (!(f_nu_p.is_hole || f_nu.is_hole)) cross_prod_2 -= (((r_nu_p - r_nu).unit())*f_nu.get_jacobian(i));
        perim_vec = perim_vec + cross_prod_2; 
        if (vi.boundary)
        {
          if (f == 0 || (f == vi.n_faces - 2))
          {
            Vector3d r_nu_i_unit = r_nu_i.unit();
            perim_vec = perim_vec + r_nu_i_unit*f_nu.get_jacobian(i) - r_nu_i_unit;
          }
        }
        //if (m_has_pair_params)
        //  lambda = m_pair_params[vi.type - 1][mesh.get_vertices()[E.to].type - 1].lambda;
        Vector3d cross_prod_3(0.0,0.0,0.0);
        if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_3 = lambda*(((r_nu - r_nu_m).unit())*f_nu.get_jacobian(i));
        con_vec = con_vec + cross_prod_3; 
        //if (m_has_pair_params)
        //  lambda = m_pair_params[vi.type - 1][mesh.get_vertices()[En.to].type - 1].lambda;
        Vector3d cross_prod_4(0.0,0.0,0.0);
        if (!(f_nu_p.is_hole || f_nu.is_hole)) cross_prod_4 = lambda*(((r_nu_p - r_nu).unit())*f_nu.get_jacobian(i));
        con_vec = con_vec - cross_prod_4; 
        if (vi.boundary)
        {
          if (f == 0 || (f == vi.n_faces - 2))
          {
            if (m_has_pair_params)
              lambda = m_pair_params[vi.type - 1][vi.type - 1].lambda;
            Vector3d r_nu_i_unit = r_nu_i.unit();
            con_vec = con_vec + lambda*(r_nu_i_unit*f_nu.get_jacobian(i) - r_nu_i_unit);
          }
        }
        pot_eng += lambda*(r_nu - r_nu_m).len();
      }
      // area term
      double area_fact = -alpha*area_term;
      pi.fx += area_fact*area_vec.x;
      pi.fy += area_fact*area_vec.y;
      pi.fz += area_fact*area_vec.z;
      if (vi.boundary)
      {
        double add_area = -2.0*area_fact*pi.A0;
        Vector3d& ang_def = vi.get_angle_def(vi.id);
        pi.fx += add_area*ang_def.x;
        pi.fy += add_area*ang_def.y;
        pi.fz += add_area*ang_def.z;
      }
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
      if (vj.attached && (m_include_boundary || !vj.boundary))  
      {
        if (m_has_part_params) 
        {
          K  = m_particle_params[vj.type-1].K;
          gamma = m_particle_params[vj.type-1].gamma;
        }
        double dA = 0.0;
        if (!vj.boundary) dA = (vj.area - pj.A0);
        else dA = (vj.area -  mesh.angle_factor(vj.id)*pj.A0);
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
            Vector3d r_nu_i = r_nu - vi.r;

            Vector3d cross_prod_1(0.0,0.0,0.0);
            if (!(f_nu.is_hole || f_nu_p.is_hole)) cross_prod_1 = cross(r_nu_p, Nvec)*f_nu.get_jacobian(i);
            if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_1 = cross_prod_1 - cross(r_nu_m, Nvec)*f_nu.get_jacobian(i);
            area_vec = area_vec + cross_prod_1; 
            //if (vi.boundary && vj.boundary)
            if (vj.boundary)
            {
               if (f == 0)
                 area_vec = area_vec - (cross(vj.r,Nvec))*f_nu.get_jacobian(i);
               else if (f == vj.n_faces - 2)
                 area_vec = area_vec + (cross(vj.r,Nvec))*f_nu.get_jacobian(i);
            }
            Vector3d cross_prod_2(0.0,0.0,0.0);
            if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_2 = ((r_nu - r_nu_m).unit())*f_nu.get_jacobian(i);
            if (!(f_nu_p.is_hole || f_nu.is_hole)) cross_prod_2 -= ((r_nu_p - r_nu).unit())*f_nu.get_jacobian(i);
            perim_vec = perim_vec + cross_prod_2; 
            if (vj.boundary)
            {
              if (f == 0 || (f == vj.n_faces - 2))
              {
                Vector3d r_nu_j = r_nu - vj.r;
                perim_vec = perim_vec + r_nu_j.unit()*f_nu.get_jacobian(i);
              }
            }
            //if (m_has_pair_params)
            //  lambda = m_pair_params[vj.type - 1][mesh.get_vertices()[E.to].type - 1].lambda;
            Vector3d cross_prod_3(0.0,0.0,0.0);
            if (!(f_nu_m.is_hole || f_nu.is_hole)) cross_prod_3 = lambda*(((r_nu - r_nu_m).unit())*f_nu.get_jacobian(i));
            con_vec = con_vec + cross_prod_3; 
            //if (m_has_pair_params)
            //  lambda = m_pair_params[vj.type - 1][mesh.get_vertices()[Ep.to].type - 1].lambda;
            Vector3d cross_prod_4(0.0,0.0,0.0);
            if (!(f_nu_p.is_hole || f_nu.is_hole)) cross_prod_4 = lambda*(((r_nu_p - r_nu).unit())*f_nu.get_jacobian(i));
            con_vec = con_vec - cross_prod_4; 
            if (vj.boundary)
            {
              if (f == 0 || (f == vj.n_faces - 2))
              {
                if (m_has_pair_params)
                  lambda = m_pair_params[vj.type - 1][vj.type - 1].lambda;
                Vector3d r_nu_j = r_nu - vj.r;
                con_vec = con_vec + lambda*(r_nu_j.unit()*f_nu.get_jacobian(i));
              }
            }
          }
        }
        // area term
        double area_fact = -alpha*area_term;
        pi.fx += area_fact*area_vec.x;
        pi.fy += area_fact*area_vec.y;
        pi.fz += area_fact*area_vec.z;
        if (vj.boundary)
        {
          double add_area = -2.0*area_fact*pj.A0;
          Vector3d& ang_def = vj.get_angle_def(vi.id);
          pi.fx += add_area*ang_def.x;
          pi.fy += add_area*ang_def.y;
          pi.fz += add_area*ang_def.z;
        }
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
