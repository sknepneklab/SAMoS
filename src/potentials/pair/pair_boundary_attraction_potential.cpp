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
 * \file pair_boundary_attraction_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Jun-2016
 * \brief Implementation of PairBoundaryAttractionPotential class
 */ 

#include "pair_boundary_attraction_potential.hpp"


//! \param dt time step sent by the integrator 
void PairBoundaryAttractionPotential::compute(double dt)
{
  double epsilon = m_epsilon;
  double rc = m_rc;
  double wc = m_wc;
  int N = m_system->size();
  Mesh& mesh = m_system->get_mesh();
  
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("boundary_attraction",0.0);
    }
  }
  
  m_potential_energy = 0.0;
  for (int v = 0; v < mesh.size(); v++)
  {
    Vertex& Vi = mesh.get_vertices()[v];
    if (Vi.boundary)
    {
      Particle& pi = m_system->get_particle(Vi.id);
      for (int j = 0; j < Vi.n_edges; j++)
      {
        Vertex& Vj = mesh.get_vertices()[Vi.neigh[j]];
        Particle& pj = m_system->get_particle(Vj.id);
        if (m_has_pair_params)
        {
          epsilon = m_pair_params[pi.get_type()-1][pj.get_type()-1].epsilon;
          rc = m_pair_params[pi.get_type()-1][pj.get_type()-1].rc;
          wc = m_pair_params[pi.get_type()-1][pj.get_type()-1].wc;
        }
        double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
        m_system->apply_periodic(dx,dy,dz);
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        if (r >= rc && r <= (rc+wc))
        {
          double cos_fac = cos(0.5*M_PI*(r-rc)/wc);
          double potential_energy = -epsilon*cos_fac*cos_fac;
          m_potential_energy += potential_energy;
          double fact = -0.5*M_PI*epsilon/wc*sin(M_PI*(r-rc)/wc)/r;
          pi.fx += fact*dx;
          pi.fy += fact*dy;
          pi.fz += fact*dz;
          pj.fx -= fact*dx;
          pj.fy -= fact*dy;
          pj.fz -= fact*dz;
          if (m_system->compute_per_particle_energy())
          {
            pi.add_pot_energy("boundary_attraction",potential_energy);
            pj.add_pot_energy("boundary_attraction",potential_energy);
          }
        }
      }
    }
  }
}