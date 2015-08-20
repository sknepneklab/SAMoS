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
 * \file log_vec_velocity.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Mar-2014
 * \brief Declaration of LogVecVelocity class
 */ 

#ifndef __LOG_VEC_VELOCITY_H__
#define __LOG_VEC_VELOCITY_H__

#include <cmath>

#include "log.hpp"

using std::sqrt;

//! LogVecVelocity class
/*! Logs magnitude of total velocity calculated as
 *  \f$ v =  \frac{1}{N}\left|\sum_i \vec v_i\right| \f$. 
 *  It is related to the Vicsek order parameters by dividing this 
 *  quantity by \f$ v_0 \f$.
 */
class LogVecVelocity : public Log
{
public:
  
  //! Construct Log object
  //! \param sys pointer to a system object
  //! \param compute pointer to a compute object
  //! \param pot Interaction handler object
  //! \param align Alignment handler object
  LogVecVelocity(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align) : Log(sys, msg, pot, align)  { }
  
  virtual ~LogVecVelocity() { }
  
   //! \return current time step
  string operator()()
  {
    double vx = 0.0, vy = 0.0, vz = 0.0;
    for (int i = 0; i < m_system->size(); i++)
    {
      Particle& p = m_system->get_particle(i);
      vx += p.vx;
      vy += p.vy;
      vz += p.vz;
    }
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    return str(format("%12.6e ") % (v/m_system->size()));
  }
  
};

typedef shared_ptr<LogVecVelocity> LogVecVelocityPtr;

#endif