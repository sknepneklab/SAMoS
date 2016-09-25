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
                       /* to add new potential: | qi::as_string[keyword["newpotential"]][phx::bind(&ExternalAlignData::type, phx::ref(external_align_data)) = qi::_1 ] */
                     )
                     >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&ExternalAlignData::params, phx::ref(external_align_data)) = qi::_1 ]
                     >> (qi::eol || qi::eoi);
  }
  
private:
  
  qi::rule<std::string::iterator, qi::space_type> external_align;  //!< Rule for parsing external alignment lines.
  
};





#endif