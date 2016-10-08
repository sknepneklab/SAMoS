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
 * \file parse_group.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Sept-2014
 * \brief Grammar for the parsing group command lines
 */ 

#ifndef __PARSE_GROUP_HPP__
#define __PARSE_GROUP_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct GroupData
{
  std::string name;     //!< contains name of the group
  std::string params;   //!< parameters that determine group, such as particle type
};

/*! This is a parser for parsing command that contain group directives.
 *
 *  For example:
 * 
 *  group small { type = 1 }
 *  group R2 { radius = 2.0 }
 * 
 * This parser will extract group name name (in this case "small" and "R2")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
 */
class group_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  group_grammar(GroupData& group_data) : group_grammar::base_type(group)
  {
    group = qi::as_string[+qi::char_("a-zA-Z_0-9.")][phx::bind(&GroupData::name, phx::ref(group_data)) = qi::_1 ]
                >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&GroupData::params, phx::ref(group_data)) = qi::_1 ]
                >> (qi::eol || qi::eoi);
  }
  
private:
  
  qi::rule<std::string::iterator, qi::space_type> group;  //!< Rule for parsing particle group lines
  
};

#endif
