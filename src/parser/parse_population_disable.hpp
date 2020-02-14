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
 * \file parse_population_disable.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-Aug-2018
 * \brief Grammar for the parsing command line that disables a population control  
 */ 

#ifndef __PARSE_DISABLE_POPULATION_HPP__
#define __PARSE_DISABLE_POPULATION_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct DisablePopulationData
{
  std::string type;      //!< population type to disable (such as nve)
  std::string params;    //!< parameters (currently not used, but might be in the future)
};

/*! This is a parser for parsing command that disable a given integrator.
 * 
 *  for example:
 * 
 *  disable_population density {  }
 * 
 * This parser disable population control of a given type (in this case "density")
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class disable_population_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  disable_population_grammar(DisablePopulationData& disable_population_data) : disable_population_grammar::base_type(disable_population)
  {
    disable_population = (
                            qi::as_string[keyword["random"]][phx::bind(&DisablePopulationData::type, phx::ref(disable_population_data)) = qi::_1 ]    /*! Disables random population */
                          | qi::as_string[keyword["density"]][phx::bind(&DisablePopulationData::type, phx::ref(disable_population_data)) = qi::_1 ]    /*! Disables density population */
                          | qi::as_string[keyword["grow"]][phx::bind(&DisablePopulationData::type, phx::ref(disable_population_data)) = qi::_1 ]    /*! Disables grow population */
                          | qi::as_string[keyword["elongate"]][phx::bind(&DisablePopulationData::type, phx::ref(disable_population_data)) = qi::_1 ]    /*! Disables elongate population */
                          | qi::as_string[keyword["cell"]][phx::bind(&DisablePopulationData::type, phx::ref(disable_population_data)) = qi::_1 ]    /*! Disables radnom population */
                          | qi::as_string[keyword["actomyosin"]][phx::bind(&DisablePopulationData::type, phx::ref(disable_population_data)) = qi::_1 ]    /*! Disables actomyosin population */
                          | qi::as_string[keyword["actomyosin_poisson"]][phx::bind(&DisablePopulationData::type, phx::ref(disable_population_data)) = qi::_1 ]    /*! Disables acotmyosin_poisson population */
                          | qi::as_string[keyword["actomyosin_molecule"]][phx::bind(&DisablePopulationData::type, phx::ref(disable_population_data)) = qi::_1 ]    /*! Disables actomyosin_molecule population */
                          | qi::as_string[keyword["actomyosin_head"]][phx::bind(&DisablePopulationData::type, phx::ref(disable_population_data)) = qi::_1 ]    /*! Disables actomyosin_head population */
                          | qi::as_string[keyword["region"]][phx::bind(&DisablePopulationData::type, phx::ref(disable_population_data)) = qi::_1 ]    /*! Disables region population */
                          /* to add new disable population: | qi::as_string[keyword["newintegrator"]][phx::bind(&DisablePopulationData::type, phx::ref(disable_population_data)) = qi::_1 ] */
                        )
                        >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&DisablePopulationData::params, phx::ref(disable_population_data)) = qi::_1 ]
                        >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> disable_population;  //!< Rule for parsing disable population lines.
  
};

#endif
