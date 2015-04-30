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
 *   (c) 2013, 2014, 2015
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file external_field_aligner.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 11-Apr-2015
 * \brief Declaration of ExternalFieldAlign class
 */ 


#include "external_field_aligner.hpp"

void ExternalFieldAlign::compute()
{
  int N = m_system->size();
  double hx = m_hx, hy = m_hy, hz = m_hz;  
  
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_has_params)
    {
      hx = m_type_params[pi.get_type()].hx;
      hy = m_type_params[pi.get_type()].hy;
      hz = m_type_params[pi.get_type()].hz;
    }
    pi.tau_x += hy*pi.nz - hz*pi.ny;
    pi.tau_y += hz*pi.nx - hx*pi.nz;
    pi.tau_z += hx*pi.ny - hy*pi.nx;
  }
}