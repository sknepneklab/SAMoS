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