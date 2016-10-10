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
 * \file constraint_plane.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Declaration of ConstraintPlane class.
 */ 

#ifndef __CONSTRAINT_PLANE_HPP__
#define __CONSTRAINT_PLANE_HPP__

#include <cmath>

#include "system.hpp"
#include "parse_parameters.hpp"
#include "constraint.hpp"

using std::sin;
using std::cos;


/*! Enforces all particles to lay on the xy plane of size Lx x Ly
 *  Also makes sure periodic boundary conditions are enforced
*/
class ConstraintPlane : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds (e.g., sphere radius)
  ConstraintPlane(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    m_msg->msg(Messenger::INFO,"Planar constraint. Setting plane dimensions to lx = "+lexical_cast<string>(m_system->get_box()->Lx)+" and ly = "+lexical_cast<string>(m_system->get_box()->Ly)+".");
    m_msg->write_config("constraint.plane.lx",lexical_cast<string>(m_system->get_box()->Lx));
    m_msg->write_config("constraint.plane.ly",lexical_cast<string>(m_system->get_box()->Ly));
    if (param.find("unlimited") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Plane constraint. For non-periodic systems particle motion will be confined by the box size.");
      m_unlimited = false;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Plane constraint. Particle motion will not be affected by the box boundaries.");
      m_unlimited = true;
    }
    if (param.find("zpos") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Plane constraint. No zpos specified. Assuming that the constraint plane is at z = 0.");
      m_zpos = 0.0;
    }
    else 
    {
      m_msg->msg(Messenger::WARNING,"Plane constraint. Setting constraint plane at z = "+param["zpos"]+".");
      m_zpos = lexical_cast<double>(param["zpos"]);
    }
    m_msg->write_config("constraint.plane.zpos",lexical_cast<string>(m_zpos));
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
  
  // Rescale constraint
  bool rescale();
   
private:
  
  bool m_unlimited;       //!< If true, ignore box boundary and low system to exapand freely
  double m_zpos;            //!< Position (along z axis) of the constraint plane
  
  
};

typedef shared_ptr<ConstraintPlane> ConstraintPlanePtr;  //!< Shared pointer to the Constraint object

#endif
