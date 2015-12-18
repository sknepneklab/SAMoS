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
 * \file log_size.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Dec-2015
 * \brief Declaration of LogSize class
 */ 

#ifndef __LOG_SIZE_H__
#define __LOG_SIZE_H__

#include "log.hpp"


//! LogSize class
/*! Logs total numner of particles in the system
 */
class LogSize : public Log
{
public:
  
  //! Construct Log object
  //! \param sys pointer to a system object
  //! \param compute pointer to a compute object
  //! \param pot Interaction handler object
  //! \param align Alignment handler object
  LogSize(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align) : Log(sys, msg, pot, align)  { }
  
  virtual ~LogSize() { }
  
   //! \return current time step
  string operator()()
  {
    return str(format("%6d ") % m_system->size());
  }
  
};

typedef shared_ptr<LogSize> LogSizePtr;

#endif