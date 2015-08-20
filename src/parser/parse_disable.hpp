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
                  qi::as_string[keyword["nve"]][phx::bind(&DisableData::type, phx::ref(disable_data)) = qi::_1 ]           /*! Disables NVE integrator */
                  | qi::as_string[keyword["nvt"]][phx::bind(&DisableData::type, phx::ref(disable_data)) = qi::_1 ]         /*! Disables NVT integrator */
                  | qi::as_string[keyword["brownian"]][phx::bind(&DisableData::type, phx::ref(disable_data)) = qi::_1 ]    /*! Disables stochastic integrator */
                  | qi::as_string[keyword["vicsek"]][phx::bind(&DisableData::type, phx::ref(disable_data)) = qi::_1 ]      /*! Disables Vicsek integrator */
                  | qi::as_string[keyword["nematic"]][phx::bind(&DisableData::type, phx::ref(disable_data)) = qi::_1 ]     /*! Disables nematic integrator */
                  /* to add new integrator: | qi::as_string[keyword["newintegrator"]][phx::bind(&DisableData::type, phx::ref(disable_data)) = qi::_1 ] */
                 )
                 >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&DisableData::params, phx::ref(disable_data)) = qi::_1 ]
                 >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> disable;  //!< Rule for parsing disable integrator lines.
  
};

#endif