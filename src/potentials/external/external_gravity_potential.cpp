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
 *   (c) 2013
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file external_gravity_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 23-Oct-2013
 * \brief Implementation of the external gravitational potential
 */ 

#include "external_gravity_potential.hpp"

/*! Apply external gravitational potential to all particles */
void compute()
{
  int N = m_system->size();
  double g = m_g;
  
  m_potential_energy = 0.0;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    if (m_has_params)
      g = m_type_params[p.type]["g"];
    m_potential_energy += g*p.z;
    p.fz -= g;
  }
}
