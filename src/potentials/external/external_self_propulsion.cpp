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
 * \file external_self_propulsion.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Dec-2015
 * \brief Implementation of the self propulsion for active particles.
 */ 

#include "external_self_propulsion.hpp"

/*! Apply self propulsion to all particles */
void ExternalSelfPropulsion::compute()
{
  int N = m_system->size();
  double alpha = m_alpha;
  
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    if (m_has_params)
      alpha = m_type_params[p.get_type()]["alpha"];
    p.fx += alpha*p.nx;
    p.fy += alpha*p.ny;
    p.fz += alpha*p.nz;
  }
}
