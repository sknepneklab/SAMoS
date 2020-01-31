/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

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
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


#include "system.hpp"
#include "parse_parameters.hpp"
#include "log.hpp"
#include "log_step.hpp"
#include "log_velocity.hpp"
#include "log_vec_velocity.hpp"
#include "log_pair_eng.hpp"
#include "log_pair_align.hpp"
#include "log_bond_eng.hpp"
#include "log_angle_eng.hpp"
#include "log_external_eng.hpp"
#include "log_area.hpp"
#include "log_perim.hpp"
#include "log_size.hpp"
#include "log_kinetic_eng.hpp"

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
