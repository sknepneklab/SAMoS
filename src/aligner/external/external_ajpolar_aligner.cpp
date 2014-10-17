/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *  
 *   Author: Silke Henkes
 * 
 *   ICSMB, Department of Physics
 *   University of Aberdeen
 * 
 *   (c) 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file external_ajpolar_aligner.cpp
 * \author Silke Henkes, silkehenkes@gmail.com
 * \date 16-Oct-2014
 * \brief Declaration of ExternalAJPolarAlign class
 */ 


#include "external_ajpolar_aligner.hpp"

void ExternalAJPolarAlign::compute()
{
  int N = m_system->size();
  double tau = m_tau;
  
  
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    // The active jamming torque is:
    // 1/tau(\theta_v-theta_i); use variant 1/tau sin(\theta_v-theta_i)
    // or in vectorial: -1/tau (n_i \times v_i). ez
    // or calculate that vector and then give it to the projection
    double tau_x = pi.ny*pi.vz - pi.nz*pi.vy;
    double tau_y = pi.nz*pi.vx - pi.nx*pi.vz;
    double tau_z = pi.nx*pi.vy - pi.ny*pi.vx;
    
    if (m_has_params)
      tau = m_type_params[pi.get_type()].tau;
    
    pi.tau_x -= tau_x/tau;
    pi.tau_y -= tau_y/tau;
    pi.tau_z -= tau_z/tau;
    
  }
}