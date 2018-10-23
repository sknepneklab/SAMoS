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
 *  disable nvt { group = all }
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
                  | qi::as_string[keyword["brownian_pos"]][phx::bind(&DisableData::type, phx::ref(disable_data)) = qi::_1 ]     /*! Disables brownian_pos integrator */
                  | qi::as_string[keyword["brownian_rod"]][phx::bind(&DisableData::type, phx::ref(disable_data)) = qi::_1 ]     /*! Disables brownian_rod integrator */
                  /* to add new integrator: | qi::as_string[keyword["newintegrator"]][phx::bind(&DisableData::type, phx::ref(disable_data)) = qi::_1 ] */
                 )
                 >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&DisableData::params, phx::ref(disable_data)) = qi::_1 ]
                 >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> disable;  //!< Rule for parsing disable integrator lines.
  
};

#endif
