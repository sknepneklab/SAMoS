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
  m_logger["step"] = std::make_shared<LogStep>(LogStep(sys, msg, pot, align));
  m_logger["velocity"] = std::make_shared<LogVelocity>(LogVelocity(sys, msg, pot, align));
  m_logger["vec_velocity"] = std::make_shared<LogVecVelocity>(LogVecVelocity(sys, msg, pot, align));
  m_logger["soft_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"soft")); 
  m_logger["soft_attractive_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"soft_attractive")); 
  m_logger["gaussian_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"gaussian"));
  m_logger["morse_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"morse")); 
  m_logger["lj_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"lj")); 
  m_logger["rod_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"rod"));
  m_logger["ljrod_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"ljrod"));
  m_logger["vp_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"vp"));
  m_logger["line_tension_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"line_tension")); 
  m_logger["boundary_bending_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"boundary_bending")); 
  m_logger["boundary_attraction_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"boundary_attraction"));
  m_logger["coulomb_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"coulomb"));
  m_logger["yukawa_energy"] = std::make_shared<LogPairEng>(LogPairEng(sys, msg, pot, align,"yukawa"));
  m_logger["polar_align"] = std::make_shared<LogPairAlign>(LogPairAlign(sys, msg, pot, align,"polar")); 
  m_logger["nematic_align"] = std::make_shared<LogPairAlign>(LogPairAlign(sys, msg, pot, align,"nematic")); 
  m_logger["harmonic_bond_energy"] = std::make_shared<LogBondEng>(LogBondEng(sys, msg, pot, align,"harmonic")); 
  m_logger["harmonic_angle_energy"] = std::make_shared<LogAngleEng>(LogAngleEng(sys, msg, pot, align,"harmonic")); 
  m_logger["cosine_angle_energy"] = std::make_shared<LogAngleEng>(LogAngleEng(sys, msg, pot, align,"cosine")); 
  m_logger["external_gravity_energy"] = std::make_shared<LogExternalEng>(LogExternalEng(sys, msg, pot, align,"gravity")); 
  m_logger["external_harmonic_energy"] = std::make_shared<LogExternalEng>(LogExternalEng(sys, msg, pot, align,"harmonic")); 
  m_logger["area"] = std::make_shared<LogArea>(LogArea(sys, msg, pot, align));
  m_logger["avg_perim"] = std::make_shared<LogPerim>(LogPerim(sys, msg, pot, align));
  m_logger["size"] = std::make_shared<LogSize>(LogSize(sys, msg, pot, align));
  m_logger["kinetic_energy"] = std::make_shared<LogKineticEng>(LogKineticEng(sys, msg, pot, align));

  // ------------------------------------------------------------------
    
  m_msg->msg(Messenger::INFO,"Added logger. Logged quantities will be sent to "+file_name+".");
  m_out.open(m_file_name.c_str()); 
  m_msg->write_config("logger."+file_name+".file_name",file_name);
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
  m_msg->write_config("logger."+file_name+".freq",lexical_cast<string>(m_freq));
  for (pairs_type::iterator it_l = params.begin(); it_l != params.end(); it_l++)
  {
    string logme = to_lower_copy((*it_l).first);
    if (logme != "freq")
    {
      if (m_logger.find(logme) != m_logger.end())
      {
        m_msg->msg(Messenger::INFO,"Adding log quantity : "+logme+".");
        m_to_log.push_back(logme);
        m_msg->add_config("logger."+file_name+".quantity",logme);
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
