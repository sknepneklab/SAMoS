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
 * \file parse_box.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Grammar for the parsing simulation box
 */ 

#ifndef __PARSE_BOX_HPP__
#define __PARSE_BOX_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct BoxData
{
  std::string type;      //!< box type (periodic or fixed)
  std::string params;    //!< box size
};

/*! This is a parser for parsing command that contain simulation box information 
 *  Structurally, it is very similar to the command parsed (\see parse_command.hpp)
 *  This parser extracts the potential line type and returns it together with 
 *  the list of parameters that will be parsed separately.
 * 
 *  for example:
 * 
 *  box periodic { lx = 10.0; ly = 10.0; lz = 1.0; }
 * 
 * This parser will extract the box type (in this case "periodic")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class box_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  box_grammar(BoxData& box_data) : box_grammar::base_type(box)
  {
    box  = (
              qi::as_string[keyword["periodic"]][phx::bind(&BoxData::type, phx::ref(box_data)) = qi::_1 ]       /*! Handles periodic box */
            | qi::as_string[keyword["fixed"]][phx::bind(&BoxData::type, phx::ref(box_data)) = qi::_1 ]          /*! Handles periodic potential */
            /* to add new box type: | qi::as_string[keyword["newboxtype"]][phx::bind(&BoxData::type, phx::ref(box_data)) = qi::_1 ] */
           )
           >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&BoxData::params, phx::ref(box_data)) = qi::_1 ]
           >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> box;  //!< Rule for parsing simulation box lines.
  
};

#endif
