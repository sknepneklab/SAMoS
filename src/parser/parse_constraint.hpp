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
 * \file parse_constraint.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Grammar for the parsing constraints
 */ 

#ifndef __PARSE_CONSTRAINT_HPP__
#define __PARSE_CONSTRAINT_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct ConstraintlData
{
  std::string type;      //!< constraint type (e.g., "sphere")
  std::string params;    //!< parameters that define constraint (e.g. R = 10)
};

/*! This is a parser for parsing command that contain imposed constraint
 * 
 *  For example:
 * 
 *  constraint sphere { R = 10.0 }
 * 
 * This parser will extract the constraint type (in this case "sphere")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
 */
class constraint_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  constraint_grammar(ConstraintlData& constraint_data) : constraint_grammar::base_type(constraint)
  {
    constraint = (
                    qi::as_string[keyword["sphere"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]      /*! Handles constraint on a sphere */
                  | qi::as_string[keyword["plane"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]       /*! Handles constraint on a plane */
                  | qi::as_string[keyword["walls"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]       /*! Handles plane walls constraints */
                  | qi::as_string[keyword["cylinder"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]    /*! Handles constraint on a cylinder */
                  | qi::as_string[keyword["peanut"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]      /*! Handles constraint on a peanut */
                  | qi::as_string[keyword["torus"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]       /*! Handles constraint on a torus */
                  | qi::as_string[keyword["ellipsoid"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]   /*! Handles constraint on an ellipsoid */
                  | qi::as_string[keyword["gyroid"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]      /*! Handles constraint on a gyroid */
                  | qi::as_string[keyword["actomyo"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]     /*! Handles constraint actomyo */
                  | qi::as_string[keyword["hourglass"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]   /*! Handles constraint on a hourglass */
                  | qi::as_string[keyword["gaussian_bump"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]   /*! Handles constraint on a Gaussian bump */
                  | qi::as_string[keyword["none"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]        /*! Handles dummy constraint */
                  | qi::as_string[keyword["tetrahedron"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ] /*! Handles constraint on a tetrahedral surface */
                  | qi::as_string[keyword["slab"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ]        /*! Handles constraint to a slab */
                  /* to add new constraint: | qi::as_string[keyword["newconstraint"]][phx::bind(&ConstraintlData::type, phx::ref(constraint_data)) = qi::_1 ] */
                 )
                 >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&ConstraintlData::params, phx::ref(constraint_data)) = qi::_1 ]
                 >> (qi::eol || qi::eoi);
  }
  
private:
  
  qi::rule<std::string::iterator, qi::space_type> constraint;  //!< Rule for parsing external potential lines.
  
};


#endif
