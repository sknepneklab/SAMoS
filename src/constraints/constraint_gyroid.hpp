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
 * \file constraint_gyroid.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Dec-2015
 * \brief Declaration of ConstraintGyroid class.
 */ 

#ifndef __CONSTRAINT_GYROID_HPP__
#define __CONSTRAINT_GYROID_HPP__

#include <cmath>

#include "system.hpp"
#include "parse_parameters.hpp"
#include "constraint.hpp"


using std::sqrt;
using std::sin;
using std::cos;

/*! Enforces all particles to be on the surface of an approximate gyroid 
 *  given by equation \f$ \cos \left(\frac{2 \pi  x}{lx}\right) \sin \left(\frac{2 \pi y}{ly}\right)+\sin \left(\frac{2 \pi  x}{lx}\right) \cos
   \left(\frac{2 \pi  z}{lz}\right)+\cos \left(\frac{2 \pi y}{ly}\right) \sin \left(\frac{2 \pi  z}{lz}\right) = 0. \f$
 * All velocities will point in the tangent direction.
*/
class ConstraintGyroid : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds 
  ConstraintGyroid(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    m_msg->write_config("constraint.gyroid.type","gyroid");
  }
  
  //! Computes normal to the surface
  void compute_normal(Particle&, double&, double&, double&); 
   
  // Computer gradient at a point
  void compute_gradient(Particle&, double&, double&, double&);
  
  // Value of the constraint
  double constraint_value(Particle&); 
    
  
};

typedef shared_ptr<ConstraintGyroid> ConstraintGyroidPtr;  //!< Shared pointer to the Constraint object

#endif
