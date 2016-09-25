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
 *   (c) 2015, 2016
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015, 2016
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */
 
/*!
 * \file external_tangent_aligner.cpp
 * \author Silke Henkes, silkehenkes@gmail.com
 * \date 16-Oct-2014
 * \brief Declaration of ExternaltangentAlign class
 */ 


#include "external_tangent_aligner.hpp"

void ExternalTangentAlign::compute()
{
  int N = m_system->size();
  double tau = m_tau;
  int Nbonds = m_system->num_bonds();
    
  for  (int i = 0; i < Nbonds; i++)
  {
    Bond& b = m_system->get_bond(i);
    Particle& pi = m_system->get_particle(b.i);
    Particle& pj = m_system->get_particle(b.j);
    double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;  // tangent along bond between i and j
    m_system->apply_periodic(dx,dy,dz);
    double l = sqrt(dx*dx + dy*dy + dz*dz);     
    dx /= l; dy /= l; dz /= l;              // normalise tangent vector
    // Compute torque between director of particle i and tangent
    double tau_x = pi.ny*dz - pi.nz*dy;
    double tau_y = pi.nz*dx - pi.nx*dz;
    double tau_z = pi.nx*dy - pi.ny*dx;
    if (m_has_params)
      tau = m_type_params[pi.get_type()].tau;
    pi.tau_x += tau_x/tau;
    pi.tau_y += tau_y/tau;
    pi.tau_z += tau_z/tau;
    // Compute torque between director of particle j and tangent
    tau_x = pj.ny*dz - pj.nz*dy;
    tau_y = pj.nz*dx - pj.nx*dz;
    tau_z = pj.nx*dy - pj.ny*dx;
    if (m_has_params)
      tau = m_type_params[pj.get_type()].tau;
    pj.tau_x += tau_x/tau;
    pj.tau_y += tau_y/tau;
    pj.tau_z += tau_z/tau;
  }
}
