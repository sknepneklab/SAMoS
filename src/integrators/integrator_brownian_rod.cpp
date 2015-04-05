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
 * \file integrator_brownian_rod.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 27-Mar-2015
 * \brief Implementation of the BrownianRod rod dynamics integrator.
 */ 

#include "integrator_brownian_rod.hpp"

/*! This function integrates equations of motion as introduced in the 
 *  "Brownian dynamics and dynamic Monte Carlo simulations of isotropic
 *   and liquid crystal phases of anisotropic colloidal particles: A comparative 
 *   study", A. Patti, A. Cuetos, Phys. Rev. E 86, 011403 (2012)
**/
void IntegratorBrownianRod::integrate()
{
  int N = m_system->get_group(m_group_name)->get_size();
  //double T = m_temp->get_val(m_system->get_run_step());
  double fd_x, fd_y, fd_z;                    // Deterministic part of the force
  vector<int> particles = m_system->get_group(m_group_name)->get_particles();
  
  // reset forces and torques
  m_system->reset_forces();
  m_system->reset_torques();
  
  // If nematic, attempt to flip directors
  if (m_nematic)
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi);
      if (m_rng->drnd() < m_tau)  // Flip direction n with probability m_tua (dt/tau, where tau is the parameter given in the input file).
      {
        p.nx = -p.nx;  p.ny = -p.ny;  p.nz = -p.nz;
      }
    }
  
  // compute forces in the current configuration
  if (m_potential)
    m_potential->compute(m_dt);
  // compute torques in the current configuration
  if (m_align)
    m_align->compute();
  // iterate over all particles 
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    double nx = p.nx, ny = p.ny, nz = p.nz;
    double a = 0.5*p.get_length();
    double b = p.get_radius();
    double denom = a*a - b*b;
    double S = 2.0/sqrt(denom)*log((a+sqrt(denom))/b);
    double D_perp = m_D0*b*((2*a*a-3*b*b)*S+2*a)/(16*M_PI*denom);
    double D_par = m_D0*b*((2*a*a-b*b)*S-2*a)/(8*M_PI*denom);
    double D_rot = 3.0*m_D0*b*((2*a*a-b*b)*S-2*a)/(16*M_PI*(a*a*a*a-b*b*b*b));
    
    // Project force onto director
    double f_dot_n = p.fx*nx + p.fy*ny + p.fz*nz;
    double fx_par = f_dot_n*nx, fy_par = f_dot_n*ny, fz_par = f_dot_n*nz;
    double fx_perp = p.fx - fx_par, fy_perp = p.fy - fy_par, fz_perp = p.fz - fz_par; 
    // compute deterministic forces
    fd_x = m_v0*p.nx + D_par*fx_par + D_perp*fx_perp; 
    fd_y = m_v0*p.ny + D_par*fy_par + D_perp*fy_perp;
    fd_z = m_v0*p.nz + D_par*fz_par + D_perp*fz_perp;
    // Update velocity
    p.vx = fd_x; 
    p.vy = fd_y;
    p.vz = fd_z;
    // Update particle position according to the eq. (1a)
    p.x += m_dt*fd_x;
    p.y += m_dt*fd_y;
    p.z += m_dt*fd_z;
    
    // Normal to the manifold
    double w1_x, w1_y, w1_z;
        
    m_constraint->compute_normal(p,w1_x,w1_y,w1_z);
    
    // normal to both director and the manifold normal, i.e., the cotangent vector
    double w2_x = ny*w1_z - nz*w1_y;
    double w2_y = nz*w1_x - nx*w1_z;
    double w2_z = nx*w1_y - ny*w1_x;
    
    double stoch = sqrt(2.0*D_rot*m_dt);
    double R1 = m_rng->gauss_rng(1.0), R2 = m_rng->gauss_rng(1.0);
    
    double dnx = p.tau_y*nz - p.tau_z*ny;
    double dny = p.tau_z*nx - p.tau_x*nz;
    double dnz = p.tau_x*ny - p.tau_y*nx;
    
    p.nx += D_rot*dnx + stoch*(R1*w1_x + R2*w2_x);
    p.ny += D_rot*dny + stoch*(R1*w1_y + R2*w2_y);
    p.nz += D_rot*dnz + stoch*(R1*w1_z + R2*w2_z);
   
    
    // Project everything back to the manifold
    m_constraint->enforce(p);
    // Update rod age   
    p.age += m_dt;
  }
}
