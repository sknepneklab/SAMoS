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
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file log_kinetic_eng.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Jul-2016
 * \brief Declaration of LogKineticEng class
 */ 

#ifndef __LOG_KINETIC_ENG_H__
#define __LOG_KINETIC_ENG_H__

#include <cmath>

#include "log.hpp"

using std::sqrt;

//! LogKineticEng class
/*! Logs kinetic energy per particle for a give time step. 
 */
class LogKineticEng : public Log
{
public:
  
  //! Construct Log object
  //! \param sys pointer to a system object
  //! \param compute pointer to a compute object
  //! \param pot Interaction handler object
  //! \param align Alignment handler object
  LogKineticEng(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align) : Log(sys, msg, pot, align)  { }
  
  virtual ~LogKineticEng() { }
  
   //! \return current time step
  string operator()()
  {
    double K = 0.0;
    for (int i = 0; i < m_system->size(); i++)
    {
      Particle& p = m_system->get_particle(i);
      double v2 = p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
      K += 0.5*v2/p.mass;
    }
    return str(format("%12.6e ") % (K/m_system->size()));
  }
  
};

typedef shared_ptr<LogKineticEng> LogKineticEngPtr;

#endif