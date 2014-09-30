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
 *   (c) 2013, 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file parse_disable.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 26-Sept-2014
 * \brief Grammar for the parsing command line that disables an integrator 
 */ 

#ifndef __PARSE_DISABLE_HPP__
#define __PARSE_DISABLE_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct DisableData
{
  std::string type;      //!< integrator type to disable (such as nve)
  std::string params;    //!< parameters (group information)
};

/*! This is a parser for parsing command that disable a given integrator.
 * 
 *  for example:
 * 
 *  disable nvt 
 * 
 * This parser disable integrator of a given type (in this case "nvt")
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class disable_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  disable_grammar(DisableData& disable_data) : disable_grammar::base_type(disable)
  {
    disable = (
                  qi::as_string[keyword["nve"]][phoenix::bind(&DisableData::type, phoenix::ref(disable_data)) = qi::_1 ]           /*! Disables NVE integrator */
                  | qi::as_string[keyword["nvt"]][phoenix::bind(&DisableData::type, phoenix::ref(disable_data)) = qi::_1 ]         /*! Disables NVT integrator */
                  | qi::as_string[keyword["brownian"]][phoenix::bind(&DisableData::type, phoenix::ref(disable_data)) = qi::_1 ]    /*! Disables stochastic integrator */
                  | qi::as_string[keyword["vicsek"]][phoenix::bind(&DisableData::type, phoenix::ref(disable_data)) = qi::_1 ]      /*! Disables Vicsek integrator */
                  | qi::as_string[keyword["nematic"]][phoenix::bind(&DisableData::type, phoenix::ref(disable_data)) = qi::_1 ]     /*! Disables nematic integrator */
                  /* to add new integrator: | qi::as_string[keyword["newintegrator"]][phoenix::bind(&DisableData::type, phoenix::ref(disable_data)) = qi::_1 ] */
                 )
                 >> qi::as_string[qi::no_skip[+qi::char_]][phoenix::bind(&DisableData::params, phoenix::ref(disable_data)) = qi::_1 ]
                 >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> disable;  //!< Rule for parsing disable integrator lines.
  
};

#endif