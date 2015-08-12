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
 * \file parse_align.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Jan-2014
 * \brief Grammar for the parsing pairwise alignment
 */ 

#ifndef __PARSE_ALIGN_HPP__
#define __PARSE_ALIGN_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct AlignData
{
  std::string type;      //!< aligner type (such as mean-field "mf")
  std::string params;    //!< global potential parameters 
};

/*! This is a parser for parsing command that contain pairwise alignment.
 *  Structurally, it is very similar to the command parsed (\see parse_command.hpp)
 *  This parser extracts the pairwise aligner line type and returns it together with 
 *  the list of parameters that will be parsed separately.
 * 
 *  for example:
 * 
 *  pair_align mf { J = 1.0; }
 * 
 * This parser will extract the pairwise aligner type (in this case "mf")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class align_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  align_grammar(AlignData& align_data) : align_grammar::base_type(align)
  {
    align = (
                qi::as_string[keyword["polar"]][phx::bind(&AlignData::type, phx::ref(align_data)) = qi::_1 ]       /*! Handles polar alignment */
              | qi::as_string[keyword["nematic"]][phx::bind(&AlignData::type, phx::ref(align_data)) = qi::_1 ]       /*! Handles nematic alignment */
              | qi::as_string[keyword["vicsek"]][phx::bind(&AlignData::type, phx::ref(align_data)) = qi::_1 ]       /*! Handles Vicsek alignment */
             /* to add new potential: | qi::as_string[keyword["newalign"]][phx::bind(&AlignData::type, phx::ref(align_data)) = qi::_1 ] */
            )
            >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&AlignData::params, phx::ref(align_data)) = qi::_1 ]
            >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> align;  //!< Rule for parsing pairwise alignment lines.
  
};

#endif