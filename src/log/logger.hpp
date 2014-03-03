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
 * \file logger.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2014
 * \brief Declaration of Logger class
 */ 

#ifndef __LOGGER_H__
#define __LOGGER_H__

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <list>
#include <boost/make_shared.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


#include "system.hpp"
#include "parse_parameters.hpp"
#include "log.hpp"
#include "log_step.hpp"
#include "log_velocity.hpp"
#include "log_pair_eng.hpp"
#include "log_pair_align.hpp"

using std::string;
using std::ofstream;
using std::map;
using std::endl;
using std::list;
using boost::format;
using boost::lexical_cast;
using boost::bad_lexical_cast;

/*! Logger class handles output of various measurements,
 *  such as average velocity, energies, etc.
 *  It is a regular text file with different 
 *  quantities printed as separate columns.
 */
class Logger
{
public:
  
  //! Constructor
  Logger(SystemPtr, MessengerPtr, PotentialPtr, AlignerPtr, const string&, pairs_type&);
  
  //! Destructor
  ~Logger() 
  {
    m_out.close();
    m_to_log.clear();
  }
  
  //! Do actual logging
  void log();
  
private:
  
  SystemPtr m_system;           //!< Pointer to the System object
  MessengerPtr m_msg;           //!< Handles system wide messages
  PotentialPtr m_potential;     //!< Handles all interactions 
  AlignerPtr m_aligner;         //!< Handles all alignment interactions 
  string m_file_name;           //!< Base file name for the output file
  pairs_type m_params;          //!< Control parameters
  int m_freq;                   //!< Frequency with which to dump data
  ofstream m_out;               //!< Output stream to which to send output
  

  map<string,LogPtr> m_logger;  //!< Contains all logger objects
  list<string> m_to_log;        //!< Holds quantities to be logged
  
  
};

typedef shared_ptr<Logger> LoggerPtr;

#endif