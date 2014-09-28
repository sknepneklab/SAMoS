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
 * \file integrator_nematic.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 01-Aug-2014
 * \brief Implementation of the nematic dynamics integrator.
 */ 

#include "integrator_nematic.hpp"

/*! This function integrates equations of motion
 *  Eqs. (1a) and (1b) of Y. Fily, et al., arXiv:1309.3714
 *  \f$ \partial_t \vec r_i = d_0 \kappa_i(t) \hat{\vec n_i} \f$ where \f$ d_0 \f$ is the 
 *  elementary displacement and \f$ \kappa_i(t) \f$ is the instantaneous value of the random 
 *  variable that can be either \f$ -1 \f$ of \f$ +1 \f$. 
 *  and
 *  \f$ \partial_t \vartheta_i = \eta_i(t) \f$, where \f$ \mu \f$ is the mobility, 
 *  \f$ \vartheta_i \f$ defines orientation of the velocity in the tangent plane,
 *  \f$ \hat{\vec n}_i = \left(\cos\vartheta_i,\sin\vartheta_i\right) \f$, and
 *  \f$ \eta_i(t) \f$ is the Gaussian white noise with zero mean and delta function 
 *  correlations, \f$ \left<\eta_i(t)\eta_j(t')\right> = 2\nu_r\delta_{ij}\delta\left(t-t'\right) \f$, with
 *  \f$ \nu_r \f$ being the rotational diffusion rate.
**/
void IntegratorNematic::integrate()
{
  int N = m_system->get_group(m_group_name)->get_size();
  vector<int> particles = m_system->get_group(m_group_name)->get_particles();
  // No need to compute potential, only compute torques in the current configuration
  if (m_align)
    m_align->compute();
  // iterate over all particles 
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    // Update particle position 
    double kappa = m_rng->drnd()-0.5;
    if (kappa != 0.0)
      kappa /= fabs(kappa);
    else
      kappa = 0.0;
    kappa *= m_d0;
    p.x += kappa*p.nx;
    p.y += kappa*p.ny;
    p.z += kappa*p.nz;
    // Project everything back to the manifold
    m_constraint->enforce(p);
    // Update angular velocity
    p.omega = m_mu*m_constraint->project_torque(p);
    //p.omega = m_dt*m_constraint->project_torque(p);
    // Change orientation of the director (in the tangent plane) according to eq. (1b)
    double dtheta = m_dt*p.omega + m_stoch_coeff*m_rng->gauss_rng(1.0);
    //double dtheta = m_dt*m_constraint->project_torque(p) + m_stoch_coeff*m_rng->gauss_rng(1.0);
    m_constraint->rotate_director(p,dtheta);
    //p.omega = dtheta*m_dt;
    p.age += m_dt;
  }
}