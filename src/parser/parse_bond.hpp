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
 * \file parse_bond.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Nov-2014
 * \brief Grammar for the parsing bond potentials
 */ 

#ifndef __PARSE_BOND_HPP__
#define __PARSE_BOND_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct BondData
{
  std::string type;      //!< bond potential type (such as harmonic or FENE)
  std::string params;    //!< global bond potential parameters 
};

/*! This is a parser for parsing command that contain bond potentials.
 *  Structurally, it is very similar to the command parsed (\see parse_command.hpp)
 *  This parser extracts the potential line type and returns it together with 
 *  the list of parameters that will be parsed separately.
 * 
 *  for example:
 * 
 *  bond harmonic { k = 1.0; l0 = 1.0 }
 * 
 * This parser will extract the bond potential type (in this case "harmonic")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class bond_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  bond_grammar(BondData& bond_data) : bond_grammar::base_type(bond)
  {
    bond = (
               qi::as_string[keyword["harmonic"]][phx::bind(&BondData::type, phx::ref(bond_data)) = qi::_1 ]       /*! Handles harmonic bonds */
             | qi::as_string[keyword["active"]][phx::bind(&BondData::type, phx::ref(bond_data)) = qi::_1 ]         /*! Handles active force bonds */
             | qi::as_string[keyword["fene"]][phx::bind(&BondData::type, phx::ref(bond_data)) = qi::_1 ]           /*! Handles FENE bonds */
            /* to add new bond potential: | qi::as_string[keyword["newpotential"]][phx::bind(&BondData::type, phx::ref(bond_data)) = qi::_1 ] */
           )
           >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&BondData::params, phx::ref(bond_data)) = qi::_1 ]
           >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> bond;  //!< Rule for parsing bond potential lines.
  
};





#endif
