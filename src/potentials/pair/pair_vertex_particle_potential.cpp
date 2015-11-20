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
    double pot_eng = 0.0;
    // We first deal with interion vertices
    int offset = 0;
    if (vi.boundary) offset = 1;    // For boundary vertices make sure to get only faces that have two neighbours
    int Nface = vi.faces.size();
    for (int f = offset; f < Nface-offset; f++)
    {
      Face& f_nu_m = mesh.get_faces()[vi.faces[prev(f,Nface)]];
      Face& f_nu   = mesh.get_faces()[vi.faces[f]];
      Face& f_nu_p = mesh.get_faces()[vi.faces[next(f,Nface)]];
      Vector3d& r_nu_m = f_nu_m.rc; 
      Vector3d& r_nu   = f_nu.rc; 
      Vector3d& r_nu_p = f_nu_p.rc; 
      double area_sum = 0.0;
      double perim_sum = 0.0;
      for (int k = 0; k < f_nu.vertices.size(); k++)
      {
        Vertex& vk = mesh.get_vertices()[f_nu.vertices[k]];
        if (!vk.boundary)
        {
          if (m_has_part_params) 
          {
            K  = m_particle_params[vk.type-1].K;
            A0 = m_particle_params[vk.type-1].A0;
            gamma = m_particle_params[vk.type-1].gamma;
          }
          area_sum += K*(vk.area - A0);
          perim_sum += gamma*vk.perim;
          pot_eng += 0.5*(K*(vk.area - A0)*(vk.area - A0)+gamma*vk.perim*vk.perim);
        }
      }
      // Area term
      Vector3d area_vec = cross(r_nu_p - r_nu_m, Nvec);
      double area_fact = 0.5*alpha*area_sum/f_nu.n_sides;
      pi.fx -= area_fact*area_vec.x;
      pi.fy -= area_fact*area_vec.y;
      pi.fz -= area_fact*area_vec.z;
      // Perimeter term
      Vector3d r1 = (r_nu - r_nu_m).unit(), r2 = (r_nu_p - r_nu).unit();
      Vector3d perim_vec = r1-r2;
      double perim_fact = alpha*perim_sum/f_nu.n_sides;
      pi.fx -= perim_fact*perim_vec.x;
      pi.fy -= perim_fact*perim_vec.y;
      pi.fz -= perim_fact*perim_vec.z;
      // contact interaction
      if (m_has_pair_params)
      {
        int e = mesh.get_edge_face()[make_pair(f_nu.id,f_nu_m.id)];
        int type_1 = vi.type;
        int type_2 = mesh.get_vertices()[mesh.get_edges()[e].other_vert(i)].type;
        lambda = m_pair_params[type_1-1][type_2-1].lambda;
      }
      double con_fact = alpha*lambda/f_nu.n_sides;
      pot_eng += lambda*(r_nu - r_nu_m).len();
      pi.fx -= con_fact*r1.x;  
      pi.fy -= con_fact*r1.y;  
      pi.fz -= con_fact*r1.z;
      if (m_has_pair_params)
      {
        int e = mesh.get_edge_face()[make_pair(f_nu_p.id,f_nu.id)];
        int type_1 = vi.type;
        int type_2 = mesh.get_vertices()[mesh.get_edges()[e].other_vert(i)].type;
        lambda = m_pair_params[type_1-1][type_2-1].lambda;
      }
      con_fact = alpha*lambda/f_nu.n_sides;
      pi.fx += con_fact*r2.x;  
      pi.fy += con_fact*r2.y;  
      pi.fz += con_fact*r2.z;
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
      pot_eng += 0.5*(K*(vi.area - A0)*(vi.area - A0)+gamma*vi.perim*vi.perim);
    }
    if (m_system->compute_per_particle_energy())
      pi.add_pot_energy("vp",pot_eng);
    m_potential_energy += pot_eng;
  }
}