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
              qi::as_string[keyword["periodic"]][phoenix::bind(&BoxData::type, phoenix::ref(box_data)) = qi::_1 ]       /*! Handles periodic box */
            | qi::as_string[keyword["fixed"]][phoenix::bind(&BoxData::type, phoenix::ref(box_data)) = qi::_1 ]          /*! Handles periodic potential */
            /* to add new box type: | qi::as_string[keyword["newboxtype"]][phoenix::bind(&BoxData::type, phoenix::ref(box_data)) = qi::_1 ] */
           )
           >> qi::as_string[qi::no_skip[+qi::char_]][phoenix::bind(&BoxData::params, phoenix::ref(box_data)) = qi::_1 ]
           >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> box;  //!< Rule for parsing simulation box lines.
  
};

#endif