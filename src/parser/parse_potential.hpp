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
 * \file parse_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 21-Oct-2013
 * \brief Grammar for the parsing pair potentials
 */ 

#ifndef __PARSE_POTENTIAL_HPP__
#define __PARSE_POTENTIAL_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct PotentialData
{
  std::string type;      //!< potential type (such as lj)
  std::string params;    //!< global potential parameters 
};

/*! This is a parser for parsing command that contain pair potentials.
 *  Structurally, it is very similar to the command parsed (\see parse_command.hpp)
 *  This parser extracts the potential line type and returns it together with 
 *  the list of parameters that will be parsed separately.
 * 
 *  for example:
 * 
 *  pair_potential lj { eps = 1.0; sigma = 1.0; r_cut = 2.5; }
 * 
 * This parser will extract the potential type (in this case "lj")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class potential_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  potential_grammar(PotentialData& potential_data) : potential_grammar::base_type(potential)
  {
    potential = (
                  qi::as_string[keyword["lj"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]         /*! Handles Lennard-Jones potential */
                  | qi::as_string[keyword["coulomb"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]  /*! Handles Coulomb potential */
                  | qi::as_string[keyword["debye"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]    /*! Handles Debye-Hueckel potential */
                  | qi::as_string[keyword["morse"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]    /*! Handles Morse potential */
                  | qi::as_string[keyword["soft"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]     /*! Handles soft-core potential */
                  | qi::as_string[keyword["gaussian"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ] /*! Handles Gaussian potential */
                  | qi::as_string[keyword["morse"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]    /*! Handles Morse potential */
                  | qi::as_string[keyword["active"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]   /*! Handles active potential */
                  | qi::as_string[keyword["rod"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]      /*! Handles soft rod potential */
                  | qi::as_string[keyword["ljrod"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]    /*! Handles Lennard-Jones rod potential */
                  | qi::as_string[keyword["soft_attractive"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]    /*! Handles soft attractive potential */
                  | qi::as_string[keyword["vp"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]       /*! Handles vertex-particle potential */
                  | qi::as_string[keyword["line_tension"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]     /*! Handles line-tension potential */
                  | qi::as_string[keyword["boundary_bending"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ] /*! Handles boundary bending potential */
                  | qi::as_string[keyword["boundary_attraction"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ] /*! Handles boundary attraction potential */
                  | qi::as_string[keyword["motor"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]    /*! Handles motor pair potential */
                  | qi::as_string[keyword["active_nematic"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]   /*! Handles active nematic potential */
                  | qi::as_string[keyword["yukawa"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]   /*! Handles Yukawa potential */
                  /* to add new potential: | qi::as_string[keyword["newpotential"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ] */
                 )
                 >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&PotentialData::params, phx::ref(potential_data)) = qi::_1 ]
                 >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> potential;  //!< Rule for parsing pair potential lines.
  
};





#endif
