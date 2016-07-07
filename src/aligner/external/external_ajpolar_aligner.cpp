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
    
      // Do normalisation here if useful
    if (m_normalise)
    {
      double vnorm=sqrt(pi.vx*pi.vx+pi.vy*pi.vy+pi.vz*pi.vz);
      if (vnorm>0){
	  tau_x /=vnorm
	  tau_y /=vnorm
	  tau_z /=vnorm
      }
      else {
	cout << "Warning, singularity in velocities ..." << endl;
    }
    
    if (m_has_params)
      tau = m_type_params[pi.get_type()].tau;
    
    pi.tau_x += tau_x/tau;
    pi.tau_y += tau_y/tau;
    pi.tau_z += tau_z/tau;
    
  }
}