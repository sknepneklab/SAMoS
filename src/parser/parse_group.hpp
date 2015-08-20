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