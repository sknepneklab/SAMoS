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
 * \file parse_log_dump.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Grammar for the parsing log and dump file command lines
 */ 

#ifndef __PARSE_LOG_DUMP_HPP__
#define __PARSE_LOG_DUMP_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct LogDumpData
{
  std::string name;     //!< contains file name for the log or dump
  std::string params;   //!< global potential parameters 
};

/*! This is a parser for parsing command that contain log or dump directives.
 *  Both cases are parsed with the same parser as both log and dump
 *  commands have identical syntax
 * 
 *  For example:
 * 
 *  log terminal { freq = 100; step; temperature; potential_energy }
 *  dump coord   { freq = 100; type = xyz; file = multi }
 * 
 * This parser will extract log or dump file name (in this case "terminal" and "coord")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
 */
class log_dump_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  log_dump_grammar(LogDumpData& log_dump_data) : log_dump_grammar::base_type(log_dump)
  {
    log_dump = qi::as_string[+qi::char_("a-zA-Z_0-9.")][phx::bind(&LogDumpData::name, phx::ref(log_dump_data)) = qi::_1 ]
                >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&LogDumpData::params, phx::ref(log_dump_data)) = qi::_1 ]
                >> (qi::eol || qi::eoi);
  }
  
private:
  
  qi::rule<std::string::iterator, qi::space_type> log_dump;  //!< Rule for parsing log and dump lines
  
};

#endif
