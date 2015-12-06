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
  double A0 = m_A0;
  double gamma = m_gamma;
  double lambda = m_lambda;
  double alpha = 1.0;  // phase in factor
  double pot_eng = 0.0;
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
    if (!vi.boundary)  // For direct intecations treat only non-boundary vertices
    {
      if (m_has_part_params) 
      {
        K  = m_particle_params[vi.type-1].K;
        A0 = m_particle_params[vi.type-1].A0;
        gamma = m_particle_params[vi.type-1].gamma;
      }
      double dA = (vi.area - A0);
      double area_term = K*dA;
      double perim_term = gamma*vi.perim;
      pot_eng = 0.5*(K*dA*dA+gamma*vi.perim*vi.perim);
      Vector3d area_vec(0.0,0.0,0.0);
      Vector3d perim_vec(0.0,0.0,0.0);
      Vector3d con_vec(0.0,0.0,0.0);
      for (int e = 0; e < vi.n_edges; e++)
      {
        Edge& Ep = mesh.get_edges()[vi.edges[prev(e,vi.n_edges)]];
        Edge& E = mesh.get_edges()[vi.edges[e]];
        Edge& En = mesh.get_edges()[vi.edges[next(e,vi.n_edges)]];
        
        Face& f_nu_m = mesh.get_faces()[Ep.face];
        Face& f_nu   = mesh.get_faces()[E.face];
        Face& f_nu_p = mesh.get_faces()[En.face];
       
        Vector3d& r_nu_m = f_nu_m.rc; 
        Vector3d& r_nu   = f_nu.rc; 
        Vector3d& r_nu_p = f_nu_p.rc; 
        area_vec = area_vec + cross(r_nu_p - r_nu_m, Nvec).scaled(0.5/f_nu.n_sides);
        perim_vec = perim_vec + ((r_nu - r_nu_m).unit()-(r_nu_p - r_nu).unit()).scaled(1.0/f_nu.n_sides);
        if (m_has_pair_params)
          lambda = m_pair_params[vi.type - 1][mesh.get_vertices()[E.to].type - 1].lambda;
        con_vec = con_vec + lambda*(r_nu - r_nu_m).unit().scaled(1.0/f_nu.n_sides);
        if (m_has_pair_params)
          lambda = m_pair_params[vi.type - 1][mesh.get_vertices()[En.to].type - 1].lambda;
        con_vec = con_vec - lambda*(r_nu_p - r_nu).unit().scaled(1.0/f_nu.n_sides);
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
    }
    // Now check neighbours
    for (int j = 0; j < vi.n_edges; j++)
    {
      Particle& pj = m_system->get_particle(vi.neigh[j]);
      Vertex& vj = mesh.get_vertices()[vi.neigh[j]];
      if (!vj.boundary)  // For direct intecations treat only non-boundary vertices
      {
        if (m_has_part_params) 
        {
          K  = m_particle_params[vj.type-1].K;
          A0 = m_particle_params[vj.type-1].A0;
          gamma = m_particle_params[vj.type-1].gamma;
        }
        double dA = (vj.area - A0);
        double area_term = K*dA;
        double perim_term = gamma*vj.perim;
        Vector3d area_vec(0.0,0.0,0.0);
        Vector3d perim_vec(0.0,0.0,0.0);
        Vector3d con_vec(0.0,0.0,0.0);
        Vector3d Nvec = Vector3d(pj.Nx, pj.Ny, pj.Nz);
        for (int e = 0; e < vj.n_edges; e++)
        {
          Edge& E = mesh.get_edges()[vj.edges[e]];
          Face& f_nu   = mesh.get_faces()[E.face];
          if (f_nu.has_vertex(i))
          {
            Edge& Ep = mesh.get_edges()[vj.edges[prev(e,vj.n_edges)]];
            Edge& En = mesh.get_edges()[vj.edges[next(e,vj.n_edges)]];
        
            Face& f_nu_m = mesh.get_faces()[Ep.face];
            Face& f_nu_p = mesh.get_faces()[En.face];
            Vector3d& r_nu_m = f_nu_m.rc; 
            Vector3d& r_nu   = f_nu.rc; 
            Vector3d& r_nu_p = f_nu_p.rc; 
            area_vec = area_vec + cross(r_nu_p - r_nu_m, Nvec).scaled(0.5/f_nu.n_sides);
            perim_vec = perim_vec + ((r_nu - r_nu_m).unit()-(r_nu_p - r_nu).unit()).scaled(1.0/f_nu.n_sides);
            if (m_has_pair_params)
              lambda = m_pair_params[vi.type - 1][mesh.get_vertices()[E.to].type - 1].lambda;
            con_vec = con_vec + lambda*(r_nu - r_nu_m).unit().scaled(1.0/f_nu.n_sides);
            if (m_has_pair_params)
              lambda = m_pair_params[vi.type - 1][mesh.get_vertices()[Ep.to].type - 1].lambda;
            con_vec = con_vec - lambda*(r_nu_p - r_nu).unit().scaled(1.0/f_nu.n_sides);
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
      }
    }
    // potential energy
    if (!vi.boundary)
    {
      if (m_has_part_params) 
      {
        K  = m_particle_params[vi.type-1].K;
        A0 = m_particle_params[vi.type-1].A0;
        gamma = m_particle_params[vi.type-1].gamma;
      }
    }
    if (m_system->compute_per_particle_energy())
      pi.add_pot_energy("vp",pot_eng);
    m_potential_energy += pot_eng;
  }
}