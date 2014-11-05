/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file bond_harmonic_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 05-Nov-2014
 * \brief Implementation of BondHarmonicPotential class
 */ 

#include "bond_harmonic_potential.hpp"

void BondHarmonicPotential::compute()
{
  int Nbonds = m_system->num_bonds();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double k = m_k;
  double l0 = m_l0;
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < Nbonds; i++)
  {
    Bond& b = m_system->get_bond(i);
    Particle& pi = m_system->get_particle(b.i);
    Particle& pj = m_system->get_particle(b.j);
    double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
    if (periodic)
    {
      if (dx > box->xhi) dx -= box->Lx;
      else if (dx < box->xlo) dx += box->Lx;
      if (dy > box->yhi) dy -= box->Ly;
      else if (dy < box->ylo) dy += box->Ly;
      if (dz > box->zhi) dz -= box->Lz;
      else if (dz < box->zlo) dz += box->Lz;
    }
    if (m_has_bond_params)
    {
      k = m_bond_params[b.type].k;
      l0 = m_bond_params[b.type].l0;
     }
      
    double r_sq = dx*dx + dy*dy + dz*dz;
    double r = sqrt(r_sq);
    double dl = r - l0;
    // Handle potential 
    double potential_energy = 0.5*k*dl*dl;
    m_potential_energy += potential_energy;
    // Handle force
    double force_factor = k*dl/r;
    pi.fx += force_factor*dx;
    pi.fy += force_factor*dy;
    pi.fz += force_factor*dz;
    // Use 3d Newton's law
    pj.fx -= force_factor*dx;
    pj.fy -= force_factor*dy;
    pj.fz -= force_factor*dz;
  }
}