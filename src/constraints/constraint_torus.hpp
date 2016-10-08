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
 * \file constraint_torus.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 15-Dec-2015
 * \brief Declaration of ConstraintTorus class.
 */ 

#ifndef __CONSTRAINT_TORUS_HPP__
#define __CONSTRAINT_TORUS_HPP__

#include <cmath>

#include "system.hpp"
#include "parse_parameters.hpp"
#include "constraint.hpp"


using std::sqrt;
using std::sin;
using std::cos;

/*! Enforces all particles to be on the surface of a torus of radius 
 *  \f$ c \f$ centred at the origin with the tube radius \f$ a \f$. All velocities will point in the tangent
 *  direction.
*/
class ConstraintTorus : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds (e.g., torus radius)
  ConstraintTorus(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    if (param.find("c") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Toroidal constraint. No torus radius set. Assuming 1.");
      m_c = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Toroidal constraint. Radius set to "+param["c"]+".");
      m_c = lexical_cast<double>(param["c"]);
    }
    m_msg->write_config("constraint.torus.c",lexical_cast<string>(m_c));
    if (param.find("a") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Toroidal constraint. No tube radius set. Assuming 0.5.");
      m_a = 0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Toroidal constraint. Tube radius set to "+param["a"]+".");
      m_a = lexical_cast<double>(param["a"]);
    }
    m_msg->write_config("constraint.torus.a",lexical_cast<string>(m_a));
    if (m_a >= m_c)
    {
      m_msg->msg(Messenger::ERROR,"Toroidal constraint. Tube radius has to be smaller than torus radius.");
      throw runtime_error("Toroidal constraint. Incompatible parameters.");
    }
  }
  
  //! Computes normal to the surface
  void compute_normal(Particle&, double&, double&, double&); 
   
  // Computer gradient at a point
  void compute_gradient(Particle&, double&, double&, double&);
  
  // Value of the constraint
  double constraint_value(Particle&); 
    
    
private:
  
  double m_a;     //!< Torus tube radius
  double m_c;     //!< Torus radius
    
};

typedef shared_ptr<ConstraintTorus> ConstraintTorusPtr;  //!< Shared pointer to the Constraint object

#endif
