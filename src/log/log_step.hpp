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
 * \file log_step.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2014
 * \brief Declaration of LogStep class
 */ 

#ifndef __LOG_STEP_H__
#define __LOG_STEP_H__

#include "log.hpp"


//! LogStep class
/*! Logs current time step
 */
class LogStep : public Log
{
public:
  
  //! Construct Log object
  //! \param sys pointer to a system object
  //! \param compute pointer to a compute object
  //! \param pot Interaction handler object
  //! \param align Alignment handler object
  LogStep(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align) : Log(sys, msg, pot, align)  { }
  
  virtual ~LogStep() { }
  
   //! \return current time step
  string operator()()
  {
    return str(format("%8d ") % m_system->get_step());
  }
  
};

typedef shared_ptr<LogStep> LogStepPtr;

#endif