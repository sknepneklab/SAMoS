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
 * \file integrator_brownian.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 01-Nov-2013
 * \brief Implementation of the Brownian dynamics integrator.
 */ 

#include "integrator_brownian.hpp"

/*! This function integrates equations of motion as introduced in the 
 *  Eqs. (1a) and (1b) of Y. Fily, et al., arXiv:1309.3714
 *  \f$ \partial_t \vec r_i = v_0 \hat{\vec n_i} + \mu \sum_{i\neq j}\vec F_{ij} \f$ and
 *  \f$ \partial_t \vartheta_i = \eta_i(t) \f$, where \f$ \mu \f$ is the mobility, 
 *  \f$ \vartheta_i \f$ defines orientation of the velocity in the tangent plane,
 *  \f$ \hat{\vec n}_i = \left(\cos\vartheta_i,\sin\vartheta_i\right) \f$, and
 *  \f$ \eta_i(t) \f$ is the Gaussian white noise with zero mean and delta function 
 *  correlations, \f$ \left<\eta_i(t)\eta_j(t')\right> = 2\nu_r\delta_{ij}\delta\left(t-t'\right) \f$, with
 *  \f$ \nu_r \f$ being the rotational diffusion rate.
**/
void IntegratorBrownian::integrate()
{
  int N = m_system->size();
  double inv_v = 1.0;
  // compute forces in the current configuration
  m_potential->compute();
  // iterate over all particles 
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    // compute unit length direction of the velocity
    double nx = p.vx, ny = p.vy, nz = p.vz;
    double len_n = sqrt(nx*nx + ny*ny + nz*nz);
    if (len_n > 0.0) inv_v = 1.0/len_n;
    nx *= inv_v;  ny *= inv_v;  nz *= inv_v;
    // Update particle position according to the eq. (1a)
    p.x += m_dt*(m_v0*nx + m_mu*p.fx);
    p.y += m_dt*(m_v0*ny + m_mu*p.fy);
    p.z += m_dt*(m_v0*nz + m_mu*p.fz);
    // Project everything back to the manifold
    m_constraint->enforce(p);
    // Change orientation of the velocity (in the tangent plane) according to eq. (1b)
    double theta = m_stoch_coeff*m_rng->gauss_rng(1.0);
    m_constraint->rotate_velocity(p,theta);
  }
}
