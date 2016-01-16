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
 *   (c) 2015, 2016
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015, 2016
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file parse_timestep.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Jan-2016
 * \brief Grammar for the parsing command that sets the global integrator time step size
 */ 

#ifndef __PARSE_TIMESTEP_HPP__
#define __PARSE_TIMESTEP_HPP__

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct TimeStepData
{
  double dt;    //!< contains integrator step size
};

/*! This is a parser for parsing command that contain number of 
 *  global step size for all integrators
 * 
 *  For example:
 * 
 *  timestep 0.001
 * 
 * This parser will return the size of the integrator step (0.001 in this case)
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class timestep_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  timestep_grammar(TimeStepData& timestep_data) : timestep_grammar::base_type(timestep)
  {
    timestep =    qi::double_[phx::bind(&TimeStepData::dt, phx::ref(timestep_data)) = qi::_1 ]
              >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> timestep;  //!< Rule for parsing timestep lines.
  
};

#endif