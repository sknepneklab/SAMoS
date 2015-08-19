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
 * \file external_harmonic_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-Aug-2015
 * \brief Implementation of the external harmonic potential
 */ 

#include "external_harmonic_potential.hpp"

/*! Apply external harmonic potential to all particles */
void ExternalHarmonicPotential::compute()
{
  int N = m_system->size();
  double k = m_k;
  double z0 = m_z0;
  
  m_potential_energy = 0.0;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    if (m_has_params)
    {
      k = m_type_params[p.get_type()]["k"];
      z0 = m_type_params[p.get_type()]["z0"];
    }
    double dz = p.z - z0;
    m_potential_energy += 0.5*k*dz*dz;
    p.fz -= k*dz;
  }
}