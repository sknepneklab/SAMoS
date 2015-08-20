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
 * \file log_external_eng.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-Aug-2015
 * \brief Declaration of LogExternalEng class
 */ 

#ifndef __LOG_EXTERNAL_ENG_H__
#define __LOG_EXTERNAL_ENG_H__


#include "log.hpp"

//! LogExternalEng class
/*! Logs different types of pair energies
 */
class LogExternalEng : public Log
{
public:
  
  //! Construct Log object
  //! \param sys pointer to a system object
  //! \param compute pointer to a compute object
  //! \param pot Interaction handler object
  //! \param align Alignment handler object
  //! \param type Type of the pair energy to log (should match the pair_potential type)
  LogExternalEng(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, const string& type) : Log(sys, msg, pot, align), m_type(type)  { }
  
  virtual ~LogExternalEng() { }
  
   //! \return current time step
  string operator()()
  {
    return str(format("%12.6e ") % m_potential->compute_external_potential_energy_of_type(m_type));
  }
  
private:

  string m_type;    //!< External energy type to log
  
};

typedef shared_ptr<LogExternalEng> LogExternalEngPtr;

#endif