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
 * \file external_ajnematic_aligner.cpp
 * \author Silke Henkes, silkehenkes@gmail.com
 * \date 16-Oct-2014
 * \brief Declaration of ExternalAJNematicAlign class
 */ 


#include "external_ajnematic_aligner.hpp"

void ExternalAJNematicAlign::compute()
{
  int N = m_system->size();
  double tau = m_tau;
  
  
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    // The active jamming torque for nematic systems is modelled as:
    // 1/tau sin(2(\theta_v-theta_i))
    // or in vectorial: -2/tau (n_i.v_i)(n_i \times v_i)
    // note: factor 2 comes from the expansion sin(2x) = 2sin(x)cos(x)
    double ni_dot_vi = pi.nx*pi.vx + pi.ny*pi.vy + pi.nz*pi.vz;
    double tau_x = pi.ny*pi.vz - pi.nz*pi.vy;
    double tau_y = pi.nz*pi.vx - pi.nx*pi.vz;
    double tau_z = pi.nx*pi.vy - pi.ny*pi.vx;
    
    if (m_has_params)
      tau = m_type_params[pi.get_type()].tau;
    
    pi.tau_x += ni_dot_vi*tau_x/tau;
    pi.tau_y += ni_dot_vi*tau_y/tau;
    pi.tau_z += ni_dot_vi*tau_z/tau;
    
  }
}