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
 * \file logger.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2014
 * \brief Definition of Logger class
 */ 

#include "logger.hpp"

/*! Constructor for a logger object (multiple different logger can be present).
 *  \param sys Pointer to the System object
 *  \param msg Global system-wide messenger 
 *  \param pot Pointer to the interaction handler
 *  \param align Pointer to the alignment handler
 *  \param file_name Name of the file to for this log
 *  \param params Various control parameters (e.g., output frequency) passed by the parser
*/ 
Logger::Logger(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, const string& file_name, pairs_type& params) :  m_system(sys), 
                                                                                                                                    m_msg(msg),
                                                                                                                                    m_potential(pot),
                                                                                                                                    m_aligner(align),
                                                                                                                                    m_file_name(file_name),
                                                                                                                                    m_params(params)
{
  
  // List of known log types
  m_logger["step"] = boost::make_shared<LogStep>(LogStep(sys, msg, pot, align));
  m_logger["velocity"] = boost::make_shared<LogVelocity>(LogVelocity(sys, msg, pot, align));
  m_logger["vec_velocity"] = boost::make_shared<LogVecVelocity>(LogVecVelocity(sys, msg, pot, align));
  m_logger["soft_energy"] = boost::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"soft")); 
  m_logger["lj_energy"] = boost::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"lj")); 
  m_logger["coulomb_energy"] = boost::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"coulomb")); 
  m_logger["polar_align"] = boost::make_shared<LogPairAlign>(LogPairAlign(sys, msg, pot, align,"polar")); 
  m_logger["nematic_align"] = boost::make_shared<LogPairAlign>(LogPairAlign(sys, msg, pot, align,"nematic")); 
  // ------------------------------------------------------------------
    
  m_msg->msg(Messenger::INFO,"Added logger. Logged quantities will be sent to "+file_name+".");
  m_out.open(m_file_name.c_str()); 
  
  // Always log step
  m_to_log.push_back("step");
  
  if (params.find("freq") == params.end())
  {
    m_msg->msg(Messenger::WARNING,"No log frequency specified. Using default of logging each 100 time steps.");
    m_freq = 100;
  }
  else
  {
    m_msg->msg(Messenger::INFO,"Logs will be produced every "+params["freq"]+" time steps.");
    m_freq = lexical_cast<int>(params["freq"]);
  }
  for (pairs_type::iterator it_l = params.begin(); it_l != params.end(); it_l++)
  {
    string logme = to_lower_copy((*it_l).first);
    if (logme != "freq")
    {
      if (m_logger.find(logme) != m_logger.end())
      {
        m_msg->msg(Messenger::INFO,"Adding log quantity : "+logme+".");
        m_to_log.push_back(logme);
      }
      else
      {
        m_msg->msg(Messenger::ERROR,"Unknown log quantity : "+logme+".");
        throw runtime_error("Unknown log quantity.");
      }
    }
  }
  ptime now = second_clock::local_time();
  m_out << "# Log file : " << m_file_name << endl;
  m_out << "# Generated on : " << now << endl;
  m_out << "# ";
  int i = 1;
  for (list<string>::iterator it_l = m_to_log.begin(); it_l != m_to_log.end(); it_l++)
    m_out << "(" << i++ << ") " << (*it_l) << " ";
  m_out << endl;
}

/*! Does actual logging
 */
void Logger::log()
{
  if ((m_system->get_step() % m_freq == 0))
  {
    for (list<string>::iterator it_l = m_to_log.begin(); it_l != m_to_log.end(); it_l++)
    {
      m_out <<  (*m_logger[*it_l])() << " ";
    }
    m_out << endl;
  }
}