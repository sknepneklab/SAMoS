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
 * \file integrator_actomyo.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 06-Nov-2014
 * \brief Implementation of the integrator for the actomyosin dynamics.
 */ 

#include "integrator_actomyo.hpp"

/*! This function integrates equations of motion for the system of
 *  actin filaments and myosin motors.
 *  Force is computed in the standard fashion. However, in addition,
 *  if \f$ \vec F_m \cdot \vec \tau > 0 \f$, where \f$ \vec F_m \f$ is the
 *  total force on myosin bead closest to the actin bead and \f$ \tau \f$
 *  if the tangent to the actin filament in the direction of the filament's
 *  head, we add a force \f$ F_a \f$ to the filament bead in the direction 
 *  of \f$ \tau \f$.
 *  \note: For simplicity, we assume that actin filament beads have type 2.
 *  This is clearly a constraint, however, makes things much simpler to 
 *  implement.
 *  \note: In the current implementation, this integrator ignores groups
 
**/
void IntegratorActomyo::integrate()
{
  int N = m_system->size();
  double T = m_temp->get_val(m_system->get_run_step()); // current temperature 
  m_stoch_coeff = sqrt(6.0*T/m_zeta*m_dt);  // I think we need 6 here (2*d) because we only project into the 2d space.
  // reset forces and torques
  m_system->reset_forces();
  m_system->reset_torques();
  
  // compute forces in the current configuration
  if (m_potential)
    m_potential->compute(m_dt);
  // iterate over all particles 
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    // If particle is actin particle (type 1) add actin force to it
    if (p.get_type() == 1)
    {
      double tx, ty, tz;
      m_system->compute_tangent(p.get_id(),tx,ty,tz);
      vector<int>& neigh = m_nlist->get_neighbours(p.get_id());
      double min_dist = 1e7;  // Some large number
      int min_j;
      bool found = false;
      for (unsigned int j = 0; j < neigh.size(); j++)
      {
        Particle& pj = m_system->get_particle(neigh[j]);
        double dx = pj.x - p.x, dy = pj.y - p.y, dz = pj.z - p.z; 
        m_system->apply_periodic(dx,dy,dz);
        double l = sqrt(dx*dx + dy*dy + dz*dz);
        if (l < min_dist && pj.get_type() == 3)  // Type 3 is active site on myosin
        {
          min_dist = l;
          min_j = j;
          found = true;
        }
      }
      if (found)
      {
        Particle& pj = m_system->get_particle(min_j);
        double dot = pj.fx*tx + pj.fy*ty + pj.fz*tz; // Check the dot product between force on myosin and tangent on actin
        if (dot > 0.0) // If it is positive, add force to the actin bead in the direction of the tangent vector
        {
          p.fx += m_f_active*tx;
          p.fy += m_f_active*ty;
          p.fz += m_f_active*tz;
        }
      }      
    }
    // Update velocity     
    p.vx = p.fx/m_zeta; 
    p.vy = p.fy/m_zeta; 
    p.vz = p.fz/m_zeta; 
    // Update particle position 
    p.x += m_dt*p.vx + m_stoch_coeff*m_rng->gauss_rng(1.0);
    p.y += m_dt*p.vy + m_stoch_coeff*m_rng->gauss_rng(1.0);
    p.z += m_dt*p.vz + m_stoch_coeff*m_rng->gauss_rng(1.0);
    // Project everything back to the manifold
    m_constraint->enforce(p);
  }
}
