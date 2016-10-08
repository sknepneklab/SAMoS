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
                  qi::as_string[keyword["gravity"]][phx::bind(&ExternalData::type, phx::ref(external_data)) = qi::_1 ]         /*! Handles gravitational potential */
                | qi::as_string[keyword["harmonic"]][phx::bind(&ExternalData::type, phx::ref(external_data)) = qi::_1 ]        /*! Handles harmonic potential */
                | qi::as_string[keyword["self_propulsion"]][phx::bind(&ExternalData::type, phx::ref(external_data)) = qi::_1 ] /*! Handles self propulsion */
                | qi::as_string[keyword["boundary_pull"]][phx::bind(&ExternalData::type, phx::ref(external_data)) = qi::_1 ]   /*! Handles boundary pull */
                  /* to add new potential: | qi::as_string[keyword["newpotential"]][phx::bind(&ExternalData::type, phx::ref(external_data)) = qi::_1 ] */
                )
                >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&ExternalData::params, phx::ref(external_data)) = qi::_1 ]
                >> (qi::eol || qi::eoi);
  }
  
private:
  
  qi::rule<std::string::iterator, qi::space_type> external;  //!< Rule for parsing external potential lines.
  
};





#endif
