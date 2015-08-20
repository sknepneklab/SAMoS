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
 * \file parse_angle.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Nov-2014
 * \brief Grammar for the parsing angle potentials
 */ 

#ifndef __PARSE_ANGLE_HPP__
#define __PARSE_ANGLE_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct AngleData
{
  std::string type;      //!< angle potential type (such as harmonic or cosine)
  std::string params;    //!< global bond potential parameters 
};

/*! This is a parser for parsing command that contain angle potentials.
 *  Structurally, it is very similar to the command parsed (\see parse_command.hpp)
 *  This parser extracts the potential line type and returns it together with 
 *  the list of parameters that will be parsed separately.
 * 
 *  for example:
 * 
 *  angle harmonic { k = 1.0; t0 = 0.0 }
 * 
 * This parser will extract the angle potential type (in this case "harmonic")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class angle_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  angle_grammar(AngleData& angle_data) : angle_grammar::base_type(angle)
  {
    angle = (
              qi::as_string[keyword["harmonic"]][phx::bind(&AngleData::type, phx::ref(angle_data)) = qi::_1 ]       /*! Handles harmonic angles */
            | qi::as_string[keyword["cosine"]][phx::bind(&AngleData::type, phx::ref(angle_data)) = qi::_1 ]         /*! Handles cosine angles */
             /* to add new angle potential: | qi::as_string[keyword["newpotential"]][phx::bind(&AngleData::type, phx::ref(angle_data)) = qi::_1 ] */
            )
            >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&AngleData::params, phx::ref(angle_data)) = qi::_1 ]
            >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> angle;  //!< Rule for parsing angle potential lines.
  
};





#endif
