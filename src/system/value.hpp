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
 * \file value.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 02-Mar-2015
 * \brief Declaration of Value class.
 */ 

#ifndef __VALUE_HPP__
#define __VALUE_HPP__

#include <string>
#include <exception>
#include <memory>

#include "messenger.hpp"
#include "parse_parameters.hpp"

using std::string;
using std::runtime_error;

using std::make_shared;
using boost::lexical_cast;
using boost::bad_lexical_cast;

//! Value class
/*! Value is the abstract class that provides interface for different forms of 
  value interpolation. For example, its children can be used to get linear 
  interpolation of the temperature between two preset values.
   
*/
class Value
{
public:
  
  //! Construct Value object
  //! \param param list of parameters that are passed to the appropriate instance
  Value(MessengerPtr msg, pairs_type& param) : m_msg(msg)
  { 
    if (param.find("min_val") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Value object. No minimum value set. Assuming 0.");
      m_min_val = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Value object. Minimum value set to "+param["min_val"]+".");
      m_min_val = lexical_cast<double>(param["min_val"]);
    }
    if (param.find("max_val") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Value object. No maximum value set. Assuming 1.0.");
      m_max_val = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Value object. Maximum value set to "+param["max_val"]+".");
      m_max_val = lexical_cast<double>(param["max_val"]);
    }
    if (param.find("steps") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Value object. No number of steps set. Assuming 1000.");
      m_steps = 1000;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Value object. Number of steps set to "+param["steps"]+".");
      m_steps = lexical_cast<int>(param["steps"]);
    }
  }
  
  virtual ~Value() { }
        
  //! Return current value
  virtual double get_val(int) = 0;
  
  
protected:
  
  MessengerPtr m_msg;         //!< Handles messages sent to output
  double m_min_val;           //!< Minimum parameter value
  double m_max_val;           //!< Maximum parameter value
  int m_steps;                //!< Number of steps between min and max values
  
};

typedef shared_ptr<Value> ValuePtr;


//! ValueConstant class
/*! ValueConstant is class returns constant value of the parameter, set to the 
 *  minimum value (min_val). All other parameters are ignored.   
*/
class ValueConstant : public Value
{
public:
  
  //! Construct ValueConstant object
  //! \param param list of parameters that are passed to the appropriate instance
  ValueConstant(MessengerPtr msg, pairs_type& param) : Value(msg, param)  
  {  
    m_msg->msg(Messenger::INFO,"Using constant value object.");
    m_msg->write_config("value.constant.val",lexical_cast<string>(m_min_val));
    m_msg->write_config("value.constant.steps",lexical_cast<string>(m_steps));
  }
        
  //! Return current value
  double get_val(int step) { return m_min_val; }
    
};

typedef shared_ptr<ValueConstant> ValueConstantPtr;

//! ValueLinear class
/*! ValueLinear uses linear interpolation to return value between min_val and
 *  max_val. If step is larger than the number of steps prescribed for the interpolation,
 *  max_val will be returned   
*/
class ValueLinear : public Value
{
public:
  
  //! Construct ValueLinear object
  //! \param param list of parameters that are passed to the appropriate instance
  ValueLinear(MessengerPtr msg, pairs_type& param) : Value(msg, param)  
  { 
    m_msg->msg(Messenger::INFO,"Using linear interpolation value object.");
    if (m_max_val < m_min_val)
    {
      m_msg->msg(Messenger::ERROR,"Value object. Minimum value has to be smaller than the maximum value.");
      throw runtime_error("Value parameters are not properly set.");
    }
    if (m_steps <= 0)
    {
      m_msg->msg(Messenger::ERROR,"Value object. Number of steps has to be positive non-zero integer.");
      throw runtime_error("Value parameters. Zero steps given.");
    }
    m_dv = (m_max_val - m_min_val)/static_cast<double>(m_steps);
    m_msg->write_config("value.linear.min_val",lexical_cast<string>(m_min_val));
    m_msg->write_config("value.linear.max_val",lexical_cast<string>(m_max_val));
    m_msg->write_config("value.linear.steps",lexical_cast<string>(m_steps));
    m_msg->write_config("value.linear.dv",lexical_cast<string>(m_dv));
  }
        
  //! Return current value
  double get_val(int step) 
  { 
    if (step > m_steps)
      return m_max_val;
    return m_min_val + static_cast<double>(step)*m_dv; 
  };
  
private:
  
  double m_dv;    //!< Value increment 
  
};

typedef shared_ptr<ValueLinear> ValueLinearPtr;


#endif
