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
 *   (c) 2013
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file parse_external.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Grammar for the parsing external potential
 */ 

#ifndef __PARSE_EXTRENAL_HPP__
#define __PARSE_EXTRENAL_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct ExternalData
{
  std::string type;      //!< potential type (such as gravity)
  std::string params;    //!< global potential parameters 
};

/*! This is a parser for parsing command that contain external potentials.
 *  Structurally, it is very similar to the pair potential parser (\see parse_potential.hpp)
 *  This parser extracts the potential line type and returns it together with 
 *  the list of parameters that will be parsed separately.
 * 
 *  for example:
 * 
 *  pair_potential gravity { g = 9.81 }
 * 
 * This parser will extract the potential type (in this case "gravity")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
 */
class external_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  external_grammar(ExternalData& external_data) : external_grammar::base_type(external)
  {
    external = (
                  qi::as_string[keyword["gravity"]][phoenix::bind(&ExternalData::type, phoenix::ref(external_data)) = qi::_1 ]       /*! Handles gravitational potential */
                  /* to add new potential: | qi::as_string[keyword["newpotential"]][phoenix::bind(&ExternalData::type, phoenix::ref(external_data)) = qi::_1 ] */
                )
                >> qi::as_string[qi::no_skip[+qi::char_]][phoenix::bind(&ExternalData::params, phoenix::ref(external_data)) = qi::_1 ]
                >> (qi::eol || qi::eoi);
  }
  
private:
  
  qi::rule<std::string::iterator, qi::space_type> external;  //!< Rule for parsing external potential lines.
  
};





#endif