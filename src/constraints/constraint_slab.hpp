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
 * \file constraint_slab.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Sep-2015
 * \brief Declaration of ConstraintSlab class.
 */ 

#ifndef __CONSTRAINT_SLAB_HPP__
#define __CONSTRAINT_SLAB_HPP__

#include <cmath>

#include "system.hpp"
#include "parse_parameters.hpp"
#include "constraint.hpp"

using std::sin;
using std::cos;


/*! Enforces all particles to lay in a slab between \f$ z_{low} \f$ and \f$ z_{hi} \f$.
 *  We assume reflective boundary at the top and the bottom
*/
class ConstraintSlab : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds (e.g., sphere radius)
  ConstraintSlab(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    if (param.find("z_lo") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Slab constraint. Lower boundary not set. Assumming 0.");
      m_z_lo = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Slab constraint. Lower boundary set to "+param["z_lo"]+".");
      m_z_lo = lexical_cast<double>(param["z_lo"]);
    }
    m_msg->write_config("constraint.slab.z_lo",lexical_cast<string>(m_z_lo));
    if (param.find("z_hi") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Slab constraint. Upper boundary not set. Assumming 1.");
      m_z_hi = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Slab constraint. Upper boundary set to "+param["z_hi"]+".");
      m_z_hi = lexical_cast<double>(param["z_hi"]);
    }
    m_msg->write_config("constraint.slab.z_hi",lexical_cast<string>(m_z_lo));
  }
  
  //! Enforce constraint
  void enforce(Particle&);
  
  //! Rotate director around normal vector to the slab (z axis)
  void rotate_director(Particle& p, double phi) { }  // Does not do anything as this is not really a surface constraint.
  
  //! Rotate velocity around normal vector to the slab (z axis)
  void rotate_velocity(Particle& p, double phi)  { } // Does not do anything as this is not really a surface constraint.
  
  //! Project torque onto normal vector to the slab (z axis) and return rotation angle change
  double project_torque(Particle& p) { return 0.0; }// Does not do anything as this is not really a surface constraint.
  
  //! Computes normal to the surface
  void compute_normal(Particle& p, double& Nx, double& Ny, double& Nz) { Nx = 0.0; Ny = 0.0; Nz = 0.0; p.Nx = Nx; p.Ny = Ny; p.Nz = Nz;  } // Does not do anything as this is not really a surface constraint.
  
  // Computer gradient at a point
  void compute_gradient(Particle& p, double& gx, double& gy, double& gz) { } // Does not do anything as this is not really a surface constraint.
  
  // Value of the constraint
  double constraint_value(Particle& p) { return 0.0; } // Does not do anything as this is not really a surface constraint.
 
private:
  
  double m_z_lo;     //!< Lower boundary of the slab
  double m_z_hi;     //!< Upper boundary of the slab

};

typedef shared_ptr<ConstraintSlab> ConstraintSlabPtr;  //!< Shared pointer to the Constraint object

#endif
