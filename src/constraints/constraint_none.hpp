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
 * \file constraint_none.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 02-Mar-2015
 * \brief Declaration of ConstraintNone class.
 */ 

#ifndef __CONSTRAINT_NONE_HPP__
#define __CONSTRAINT_NONE_HPP__

#include <cmath>

#include "system.hpp"
#include "parse_parameters.hpp"
#include "constraint.hpp"


/*! A dummy constraint that does nothing, but has to be there for
 *  technical reasons
*/
class ConstraintNone : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds (e.g., sphere radius)
  ConstraintNone(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    m_msg->msg(Messenger::INFO,"Constraint none.");
  }
  
  //! Enforce constraint
  void enforce(Particle& p) 
  { 
    if (m_system->get_periodic())
      m_system->enforce_periodic(p);
  }
  
  //! Rotate director around normal vector to the plane (z axis)
  void rotate_director(Particle& p, double phi) { }
  
  //! Rotate velocity around normal vector to the plane (z axis)
  void rotate_velocity(Particle& p, double phi) { }
  
  //! Project torque onto normal vector to the plane (z axis) and return rotation angle change
  double project_torque(Particle& p) { return 0.0; }
  
  //! Computes normal to the surface
  void compute_normal(Particle& p, double& Nx, double& Ny, double& Nz) { p.Nx = 0.0; p.Ny = 0.0; p.Nz = 0.0; }
  
  // Computer gradient at a point
  void compute_gradient(Particle& p, double& gx, double& gy, double& gz) { }
  
  // Value of the constraint
  double constraint_value(Particle& p) { return 0.0; }
    
  
};

typedef shared_ptr<ConstraintNone> ConstraintNonePtr;  //!< Shared pointer to the Constraint object

#endif
