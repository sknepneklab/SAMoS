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
 * \file parse_integrator.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Grammar for the parsing integrator
 */ 

#ifndef __PARSE_INTEGRATOR_HPP__
#define __PARSE_INTEGRATOR_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct IntegratorData
{
  std::string type;      //!< integrator type (such as nve)
  std::string params;    //!< integrator parameters
};

/*! This is a parser for parsing command that contain controls the integrator.
 *  Structurally, it is very similar to the potential parser (\see parse_potential.hpp)
 *  This parser extracts the integrator type and returns it together with 
 *  the list of parameters that will be parsed separately.
 * 
 *  for example:
 * 
 *  integrator nvt { dt = 0.005; Tstart = 1.0; Tend = 0.5; tau = 0.5;  }
 * 
 * This parser will integrator type (in this case "nvt")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class integrator_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  integrator_grammar(IntegratorData& integrator_data) : integrator_grammar::base_type(integrator)
  {
    integrator = ( 
                  qi::as_string[keyword["nve"]][phx::bind(&IntegratorData::type, phx::ref(integrator_data)) = qi::_1 ]              /*! Handles NVE integrator */
                  | qi::as_string[keyword["nvt"]][phx::bind(&IntegratorData::type, phx::ref(integrator_data)) = qi::_1 ]            /*! Handles NVT integrator */
                  | qi::as_string[keyword["brownian"]][phx::bind(&IntegratorData::type, phx::ref(integrator_data)) = qi::_1 ]       /*! Handles stochastic integrator */
                  | qi::as_string[keyword["vicsek"]][phx::bind(&IntegratorData::type, phx::ref(integrator_data)) = qi::_1 ]         /*! Handles Vicsek integrator */
                  | qi::as_string[keyword["nematic"]][phx::bind(&IntegratorData::type, phx::ref(integrator_data)) = qi::_1 ]        /*! Handles nematic integrator */
                  | qi::as_string[keyword["actomyo"]][phx::bind(&IntegratorData::type, phx::ref(integrator_data)) = qi::_1 ]        /*! Handles actomyo integrator */
                  | qi::as_string[keyword["brownian_rod"]][phx::bind(&IntegratorData::type, phx::ref(integrator_data)) = qi::_1 ]   /*! Handles stochastic integrator for rods */
                  | qi::as_string[keyword["brownian_pos"]][phx::bind(&IntegratorData::type, phx::ref(integrator_data)) = qi::_1 ]   /*! Handles stochastic integrator for particle position */
                  | qi::as_string[keyword["brownian_align"]][phx::bind(&IntegratorData::type, phx::ref(integrator_data)) = qi::_1 ] /*! Handles stochastic integrator for particle alignment */
                  | qi::as_string[keyword["langevin"]][phx::bind(&IntegratorData::type, phx::ref(integrator_data)) = qi::_1 ]       /*! Handles Langevin stochastic integrator */
                  /* to add new integrator: | qi::as_string[keyword["newintegrator"]][phx::bind(&IntegratorData::type, phx::ref(integrator_data)) = qi::_1 ] */
                 )
                 >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&IntegratorData::params, phx::ref(integrator_data)) = qi::_1 ]
                 >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> integrator;  //!< Rule for parsing integrator lines.
  
};

#endif