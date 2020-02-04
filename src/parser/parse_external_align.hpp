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
 * \file parse_external_align.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Jan-2014
 * \brief Grammar for the parsing external alignment
 */ 

#ifndef __PARSE_EXTRENAL_ALIGN_HPP__
#define __PARSE_EXTRENAL_ALIGN_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct ExternalAlignData
{
  std::string type;      //!< external aligner type (such as gravity)
  std::string params;    //!< global parameters 
};

/*! This is a parser for parsing command that contain external alignment.
 *  Structurally, it is very similar to the pair potential parser (\see parse_potential.hpp)
 *  This parser extracts the external alignment line type and returns it together with 
 *  the list of parameters that will be parsed separately.
 * 
 *  for example:
 * 
 *  external_align gravity { g = 9.81 }
 * 
 * This parser will extract the aligner type (in this case "gravity")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
 */
class external_align_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  external_align_grammar(ExternalAlignData& external_align_data) : external_align_grammar::base_type(external_align)
  {
    external_align = (
                       qi::as_string[keyword["gravity"]][phx::bind(&ExternalAlignData::type, phx::ref(external_align_data)) = qi::_1 ]  /*! Handles gravitational alignment */
                     | qi::as_string[keyword["ajpolar"]][phx::bind(&ExternalAlignData::type, phx::ref(external_align_data)) = qi::_1 ]       /*! Handles active jamming polar alignment */
                     | qi::as_string[keyword["ajnematic"]][phx::bind(&ExternalAlignData::type, phx::ref(external_align_data)) = qi::_1 ]       /*! Handles active jamming nematic alignment */
                     | qi::as_string[keyword["field"]][phx::bind(&ExternalAlignData::type, phx::ref(external_align_data)) = qi::_1 ]    /*! Handles alignment to external vector field */
                     | qi::as_string[keyword["cell_shape"]][phx::bind(&ExternalAlignData::type, phx::ref(external_align_data)) = qi::_1 ]    /*! Handles alignment to cell shape */
                     | qi::as_string[keyword["tangent"]][phx::bind(&ExternalAlignData::type, phx::ref(external_align_data)) = qi::_1 ]    /*! Handles alignment polymer tangent */
                     | qi::as_string[keyword["kenotaxis"]][phx::bind(&ExternalAlignData::type, phx::ref(external_align_data)) = qi::_1 ]    /*! Handles kenotaxis alignment  */
                     | qi::as_string[keyword["radial"]][phx::bind(&ExternalAlignData::type, phx::ref(external_align_data)) = qi::_1 ]    /*! Handles radial alignment  */
                     | qi::as_string[keyword["piv"]][phx::bind(&ExternalAlignData::type, phx::ref(external_align_data)) = qi::_1 ]    /*! Handles radial alignment  */
                       /* to add new potential: | qi::as_string[keyword["newpotential"]][phx::bind(&ExternalAlignData::type, phx::ref(external_align_data)) = qi::_1 ] */
                     )
                     >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&ExternalAlignData::params, phx::ref(external_align_data)) = qi::_1 ]
                     >> (qi::eol || qi::eoi);
  }
  
private:
  
  qi::rule<std::string::iterator, qi::space_type> external_align;  //!< Rule for parsing external alignment lines.
  
};





#endif
