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
 * \file pair_line_tension_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 07-Dec-2015
 * \brief Implementation of PairLineTensionPotential class
 */ 

#include "pair_line_tension_potential.hpp"


//! \param dt time step sent by the integrator 
void PairLineTensionPotential::compute(double dt)
{
  double lambda = m_lambda;
  int N = m_system->size();
  Mesh& mesh = m_system->get_mesh();
  
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("line_tension",0.0);
    }
  }
  
  m_potential_energy = 0.0;
  for (int e = 0; e < mesh.nedges(); e++)
  {
    Edge& E = mesh.get_edges()[e];
    if (E.boundary)
    {
      Vertex& Vi = mesh.get_vertices()[E.from];
      Vertex& Vj = mesh.get_vertices()[E.to];
      Particle& pi = m_system->get_particle(E.from);
      Particle& pj = m_system->get_particle(E.to);
      if (m_has_pair_params)
        lambda = m_pair_params[pi.get_type() - 1][pj.get_type() - 1].lambda;
      Vector3d rij = Vj.r - Vi.r;
      double l = rij.len();
      Vector3d urij = rij.unit();
      double fact = lambda*l;
      pi.fx += fact*urij.x;
      pi.fy += fact*urij.y;
      pi.fz += fact*urij.z;
      pj.fx -= fact*urij.x;
      pj.fy -= fact*urij.y;
      pj.fz -= fact*urij.z;
      double pot_eng = 0.5*lambda*l*l;
      m_potential_energy += pot_eng;
      if (m_system->compute_per_particle_energy())
      {
        pi.add_pot_energy("line_tension",pot_eng);
        pj.add_pot_energy("line_tension",pot_eng);
      }
    }
  }
}