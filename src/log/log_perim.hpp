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
 * \file log_perim.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Dec-2015
 * \brief Declaration of LogPerim class
 */ 

#ifndef __LOG_PERIM_H__
#define __LOG_PERIM_H__

#include "log.hpp"


//! LogPerim class
/*! Logs average perimeter for a cell system
 */
class LogPerim : public Log
{
public:
  
  //! Construct Log object
  //! \param sys pointer to a system object
  //! \param compute pointer to a compute object
  //! \param pot Interaction handler object
  //! \param align Alignment handler object
  LogPerim(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align) : Log(sys, msg, pot, align)  { }
  
  virtual ~LogPerim() { }
  
   //! \return current time step
  string operator()()
  {
    return str(format("%10.8f ") % m_system->compute_average_perimeter());
  }
  
};

typedef shared_ptr<LogPerim> LogPerimPtr;

#endif