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
 * \file constraint_ellipsoid.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Dec-2015
 * \brief Declaration of ConstraintEllipsoid class.
 */ 

#ifndef __CONSTRAINT_ELLIPSOID_HPP__
#define __CONSTRAINT_ELLIPSOID_HPP__

#include <cmath>

#include "system.hpp"
#include "parse_parameters.hpp"
#include "constraint.hpp"


using std::sqrt;
using std::sin;
using std::cos;

/*! Enforces all particles to be on the surface of an ellipsoid with parameters 
 *  \f$ a \f$, \f$ b \f$ and \f$ c \f$ All velocities will point in the tangent
 *  direction.
*/
class ConstraintEllipsoid : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds 
  ConstraintEllipsoid(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    if (param.find("a") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Ellipsoid constraint. No parameter a. Assuming 1.");
      m_a = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Ellipsoid constraint. Parameter a set to "+param["a"]+".");
      m_a = lexical_cast<double>(param["a"]);
    }
    m_msg->write_config("constraint.ellipsoid.a",lexical_cast<string>(m_a));
    if (param.find("b") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Ellipsoid constraint. No parameter b. Assuming 1.");
      m_b = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Ellipsoid constraint. Parameter b set to "+param["b"]+".");
      m_b = lexical_cast<double>(param["b"]);
    }
    m_msg->write_config("constraint.ellipsoid.b",lexical_cast<string>(m_b));
    if (param.find("c") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Ellipsoid constraint. No parameter c. Assuming 1.");
      m_c = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Ellipsoid constraint. Parameter c set to "+param["c"]+".");
      m_c = lexical_cast<double>(param["c"]);
    }
    m_msg->write_config("constraint.ellipsoid.c",lexical_cast<string>(m_c));
  }
  
  //! Computes normal to the surface
  void compute_normal(Particle&, double&, double&, double&); 
   
  // Computer gradient at a point
  void compute_gradient(Particle&, double&, double&, double&);
  
  // Value of the constraint
  double constraint_value(Particle&); 
    
private:
  
  double m_a;     //!< Ellipsoid parameter a
  double m_b;     //!< Ellipsoid parameter b
  double m_c;     //!< Ellipsoid parameter c
    
};

typedef shared_ptr<ConstraintEllipsoid> ConstraintEllipsoidPtr;  //!< Shared pointer to the Constraint object

#endif
