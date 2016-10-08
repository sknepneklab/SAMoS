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
 * \file parse_input.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Grammar for the parsing input file line
 */ 

#ifndef __PARSE_INPUT_HPP__
#define __PARSE_INPUT_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct InputData
{
  std::string name;    //!< contains input file name
};

/*! This is a parser for parsing command that contain input file
 *  generator seed.
 * 
 *  For example:
 * 
 *  input test.dat
 * 
 * This parser will return the file name ('test.dat' in this case)
 * 
 * \note We use the same parser to parse internal messages command such as in 
 * messages terminal 
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class input_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  input_grammar(InputData& input_data) : input_grammar::base_type(input)
  {
    input = qi::as_string[+qi::char_][phx::bind(&InputData::name, phx::ref(input_data)) = qi::_1 ]
            >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> input;  //!< Rule for input filename line
  
};

#endif
