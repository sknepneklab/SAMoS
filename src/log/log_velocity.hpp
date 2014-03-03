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
 * \file log_velocity.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2014
 * \brief Declaration of LogVelocity class
 */ 

#ifndef __LOG_VELOCITY_H__
#define __LOG_VELOCITY_H__

#include <cmath>

#include "log.hpp"

using std::sqrt;

//! LogVelocity class
/*! Logs magnitude of total velocity 
 */
class LogVelocity : public Log
{
public:
  
  //! Construct Log object
  //! \param sys pointer to a system object
  //! \param compute pointer to a compute object
  //! \param pot Interaction handler object
  //! \param align Alignment handler object
  LogVelocity(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align) : Log(sys, msg, pot, align)  { }
  
  virtual ~LogVelocity() { }
  
   //! \return current time step
  string operator()()
  {
    double avr_v = 0.0;
    for (int i = 0; i < m_system->size(); i++)
    {
      Particle& p = m_system->get_particle(i);
      avr_v += sqrt(p.vx*p.vx + p.vy*p.vy + p.vz*p.vz);
    }
    return str(format("%12.6e ") % (avr_v/m_system->size()));
  }
  
};

typedef shared_ptr<LogVelocity> LogVelocityPtr;

#endif