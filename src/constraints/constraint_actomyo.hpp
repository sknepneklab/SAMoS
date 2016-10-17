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
 * \file constraint_actomyo.hpp
 * \author Amit Das, dosamit@gmail.com
 * \date 2--Nov-2014
 * \brief Declaration of ConstraintActomyo class.
 */ 

#ifndef __CONSTRAINT_ACTOMYO_HPP__
#define __CONSTRAINT_ACTOMYO_HPP__

#include <cmath>

#include "system.hpp"
#include "parse_parameters.hpp"
#include "constraint.hpp"

using std::sin;
using std::cos;


/*! Enforces all particles of a actin to lay on a plane parallel to the xy direction of size Lx x Ly and
 *   particles of myosin to move at another plane of same size above the former containing the actins
 *  Also makes sure periodic boundary conditions are enforced
*/
class ConstraintActomyo : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds (e.g., sphere radius)
  ConstraintActomyo(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    m_lx = m_system->get_box()->Lx;
    m_ly = m_system->get_box()->Ly;
    m_lz = m_system->get_box()->Lz;
    m_msg->msg(Messenger::INFO,"Actomyosin constraint. Setting box dimensions to lx = "+lexical_cast<string>(m_lx)+", ly = "+lexical_cast<string>(m_ly)+", and lz ="+lexical_cast<string>(m_lz)+".");
    m_msg->write_config("constraint.actomyo.lx",lexical_cast<string>(m_lx));
    m_msg->write_config("constraint.actomyo.ly",lexical_cast<string>(m_ly));
    m_msg->write_config("constraint.actomyo.lz",lexical_cast<string>(m_lz));
  }
  
  //! Enforce constraint
  void enforce(Particle&);
  
  //! Rotate director around normal vector to the plane (z axis)
  void rotate_director(Particle&, double);
  
  //! Rotate velocity around normal vector to the plane (z axis)
  void rotate_velocity(Particle&, double);
  
  //! Project torque onto normal vector to the plane (z axis) and return rotation angle change
  double project_torque(Particle&);
  
  //! Computes normal to the surface
  void compute_normal(Particle& p, double& Nx, double& Ny, double& Nz) { Nx = 0.0; Ny = 0.0; Nz = 1.0; p.Nx = Nx; p.Ny = Ny; p.Nz = Nz; }
  
  // Computer gradient at a point
  void compute_gradient(Particle& p, double& gx, double& gy, double& gz) { }
  
  // Value of the constraint
  double constraint_value(Particle& p) { return 0.0; }
    
private:
  
  double m_lx;     //!< box size in x direction 
  double m_ly;     //!< box size in y direction 
  double m_lz;     //!< box size in z direction 
  
};

typedef shared_ptr<ConstraintActomyo> ConstraintActomyoPtr;  //!< Shared pointer to the Constraint object

#endif
