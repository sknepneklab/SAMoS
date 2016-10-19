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
 * \file messenger.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Oct-2013
 * \brief Definition of Messenger class
 */ 

#include "messenger.hpp"

//! Construct Messenger object
//! \param name file name to which to send messages or "terminal" to send it to terminal
Messenger::Messenger(const string& name) : m_file_name(name), m_has_config(false) //, m_out( (to_lower_copy(name) == "terminal")  ? cout : *(new ofstream(m_file_name.c_str())) )
{ 
  if (to_lower_copy(name) != "terminal")
  {
    m_out.open(m_file_name.c_str());
    m_to_terminal = false;
  }
  else m_to_terminal = true;
  m_msg_type[Messenger::INFO] = "INFO: ";
  m_msg_type[Messenger::WARNING] = "WARNING: ";
  m_msg_type[Messenger::ERROR] = "ERROR: ";
}

//! Destructor
Messenger::~Messenger() 
{
  if (!m_to_terminal)
    m_out.close();
  if (m_has_config)
  {
    if (m_config_type == "json")
      pt::write_json(m_config_file+".json",m_config);
    else if (m_config_type == "info")
      pt::write_info(m_config_file+".info",m_config);
    else 
      pt::write_xml(m_config_file+".xml",m_config);
  }
}

//! output message
//! Message format: TYPE: DATE TIME: MESSAGE
//! \param type Type of the message (INFO, WARNING, or ERROR)
//! \param m The message to print
void Messenger::msg(const Messenger::MSG_TYPE& type,const string& m)
{
  ptime now = second_clock::local_time();
  string val = m_msg_type.find(type)->second; // Workaround the fact that std::map::operator[] cannot pass constant arguments
  if (m_to_terminal)
    cout << val << now << ": " << m << endl;
  else
    m_out << val << now << ": " << m << endl;
}
