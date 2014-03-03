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
 * \file log.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2014
 * \brief Declaration of Log class
 */ 

#ifndef __LOG_H__
#define __LOG_H__

#include <string>

#include "system.hpp"
#include "potential.hpp"
#include "aligner.hpp"

using std::string;


//! Log class
/*! This is the base abstract class for logs of different quantities.
 *  It provides an abstract virtual operator() that returns quantity 
 *  converted to string. Children will implement it.
 */
class Log
{
public:
  
  //! Construct Log object
  //! \param sys pointer to a system object
  //! \param compute pointer to a compute object
  //! \param pot Interaction handler object
  //! \param align Alignment handler object
  Log(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align) : m_system(sys), m_msg(msg), m_potential(pot), m_aligner(align) {}
  
  virtual ~Log() { }
  
  //! Makes Log a functor
  //! \return logged quantity
  virtual string operator()() = 0;
  
protected:
  
  SystemPtr m_system;        //!< Pointer to System object
  MessengerPtr m_msg;        //!< Handles system wide messages
  PotentialPtr m_potential;  //!< Pointer to the interaction handler
  AlignerPtr m_aligner;      //!< Pointer to the alignment handler
  
};

typedef shared_ptr<Log> LogPtr;

#endif