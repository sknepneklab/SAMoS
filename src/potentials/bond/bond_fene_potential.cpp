/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
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
 * \file bond_fene_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Sept-2015
 * \brief Implementation of BondFenePotential class
 */ 

#include "bond_fene_potential.hpp"

void BondFenePotential::compute()
{
  int Nbonds = m_system->num_bonds();
  BoxPtr box = m_system->get_box();
  double k = m_k;
  double r0 = m_r0;
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < Nbonds; i++)
  {
    Bond& b = m_system->get_bond(i);
    Particle& pi = m_system->get_particle(b.i);
    Particle& pj = m_system->get_particle(b.j);
    double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
    m_system->apply_periodic(dx,dy,dz);
    if (m_has_bond_params)
    {
      k = m_bond_params[b.type-1].k;
      r0 = m_bond_params[b.type-1].r0;
    }
      
    double r_sq = dx*dx + dy*dy + dz*dz;
    // Handle potential 
    double r0_sq = r0*r0;
    double fact = 1.0-r_sq/r0_sq;
    double potential_energy = -0.5*k*r0_sq*log(fact);
    m_potential_energy += potential_energy;
    // Handle force
    double force_factor = -k/fact;
    pi.fx += force_factor*dx;
    pi.fy += force_factor*dy;
    pi.fz += force_factor*dz;
    // Use 3d Newton's law
    pj.fx -= force_factor*dx;
    pj.fy -= force_factor*dy;
    pj.fz -= force_factor*dz;
  }
}