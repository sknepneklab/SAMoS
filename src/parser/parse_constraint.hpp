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
  constraint_grammar(ConstraintlData& constraint_data) : potential_grammar::base_type(constraint)
  {
    constraint = (
                   qi::as_string[keyword["sphere"]][phoenix::bind(&ConstraintlData::type, phoenix::ref(constraint_data)) = qi::_1 ]       /*! Handles constraint on a sphere */
                   /* to add new constraint: | qi::as_string[keyword["newconstraint"]][phoenix::bind(&ConstraintlData::type, phoenix::ref(constraint_data)) = qi::_1 ] */
                 )
                 >> qi::as_string[qi::no_skip[+qi::char_]][phoenix::bind(&ConstraintlData::params, phoenix::ref(constraint_data)) = qi::_1 ]
                 >> (qi::eol || qi::eoi);
  }
  
private:
  
  qi::rule<std::string::iterator, qi::space_type> constraint;  //!< Rule for parsing external potential lines.
  
};


#endif