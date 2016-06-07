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
  double l0 = m_l0;
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
      Particle& pi = m_system->get_particle(E.from);
      Particle& pj = m_system->get_particle(E.to);
      double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
      m_system->apply_periodic(dx,dy,dz);
      if (m_has_pair_params)
      {
        lambda = m_pair_params[pi.get_type() - 1][pj.get_type() - 1].lambda;
        l0 = m_pair_params[pi.get_type() - 1][pj.get_type() - 1].l0;
      }
      double r_sq = dx*dx + dy*dy + dz*dz;
      double r = sqrt(r_sq);
      double dl = r - l0;
      double pot_eng = 0.5*lambda*dl*dl;
      m_potential_energy += pot_eng;
      // Handle force
      double force_factor = lambda*dl/r;
      pi.fx += force_factor*dx;
      pi.fy += force_factor*dy;
      pi.fz += force_factor*dz;
      // Use 3d Newton's law
      pj.fx -= force_factor*dx;
      pj.fy -= force_factor*dy;
      pj.fz -= force_factor*dz;
      if (m_system->compute_per_particle_energy())
      {
        pi.add_pot_energy("line_tension",pot_eng);
        pj.add_pot_energy("line_tension",pot_eng);
      }
    }
  }
}