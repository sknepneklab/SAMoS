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
 * \file constraint.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 24-Oct-2013
 * \brief Declaration of Constraint class.
 */ 

#ifndef __CONSTRAINT_HPP__
#define __CONSTRAINT_HPP__

#include "defaults.hpp"

#include "system.hpp"

#include "parse_parameters.hpp"


/*! Base class for a generic constraint
 *  \note This is not a constraint (such as fixed volume, or bond lengths)
 *  in the sense as it is usually defined in standard texts on molecular dynamics, 
 *  but it represents an enforcement for the particles to sit on the surface 
 *  of a manifold (e.g., sphere) and velocities to be tangent to the manifold
 *  at every given point.
*/
class Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds (e.g., sphere radius)
  Constraint(SystemPtr sys, MessengerPtr msg, pairs_type& param) : m_system(sys),
                                                                   m_msg(msg), 
                                                                   m_rescale(1.0),
                                                                   m_rescale_steps(1000),
                                                                   m_rescale_freq(10),
                                                                   m_group("all") 
  { 
    if (param.find("maxiter") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Constraint. Maximum number of iterations has not been set. Assuming 100");
      m_max_iter = 100;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Constraint. Maximum number of iterations set to "+param["maxiter"]+".");
      m_max_iter = lexical_cast<int>(param["maxiter"]);
    }
    m_msg->write_config("constraint.max_iter",lexical_cast<string>(m_max_iter));
    if (param.find("tol") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Constraint. Tolerance has not been set. Assuming 1e-6.");
      m_tol = 1e-6;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Constraint. Tolerance set to "+param["tol"]+".");
      m_tol = lexical_cast<int>(param["tol"]);
    }
    m_msg->write_config("constraint.tol",lexical_cast<string>(m_tol));
    if (param.find("rescale") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Constraint. Rescaling constraint by factor "+param["rescale"]+".");
      m_rescale = lexical_cast<double>(param["rescale"]);
      m_msg->write_config("constraint.rescale","true");
    }
    if (param.find("rescale_steps") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Constraint. Rescaling of the constraint done over "+param["rescale_steps"]+" time steps.");
      m_rescale_steps = lexical_cast<double>(param["rescale_steps"]);
      m_msg->write_config("constraint.rescale_steps",lexical_cast<string>(m_rescale_steps));
    }
    if (param.find("rescale_freq") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Constraint. Rescaling constraint every "+param["rescale_freq"]+" time steps.");
      m_rescale_freq = lexical_cast<double>(param["rescale_freq"]);
      m_msg->write_config("constraint.rescale_frequency",lexical_cast<string>(m_rescale_freq));
    }
    m_scale = pow(m_rescale,static_cast<double>(m_rescale_freq)/static_cast<double>(m_rescale_steps));
    if (m_scale != 1.0)
    {
      m_msg->msg(Messenger::INFO,"Constraint. Rescaling constraint with scale factor"+lexical_cast<string>(m_scale)+" per time steps.");
      m_msg->write_config("constraint.scale",lexical_cast<string>(m_scale));
    }
    if (param.find("group") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Constraint. No particle groups set. Assuming \"all\".");
      m_group = "all";
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Constraint. Applying constraint go group "+param["group"]+".");
      m_group = param["group"];
    }
    m_msg->write_config("constraint.group",m_group);
  }
  
  //! Enforce constraint
  virtual void enforce(Particle&);
  
  //! Rotate director around normal vector to the surface
  virtual void rotate_director(Particle&, double);
  
  //! Rotate velocity around normal vector to the surface
  virtual void rotate_velocity(Particle&, double);
  
  //! Project torque onto normal vector and return rotation angle change
  virtual double project_torque(Particle&);
  
  //! Computes normal to the surface
  virtual void compute_normal(Particle&, double&, double&, double&) = 0;
  
  // Computer gradient at a point
  virtual void compute_gradient(Particle&, double&, double&, double&) = 0;
  
  // Value of the constraint
  virtual double constraint_value(Particle&) = 0;
  
  // Rescale constraint
  virtual bool rescale() { return false;}
  
  //! Return the constraint group
  string get_group() { return m_group; }
  
protected:
  
  SystemPtr  m_system;              //!< Pointer to the system object
  MessengerPtr m_msg;               //!< Handles internal messages
  int m_max_iter;                   //!< Maximum number of iterations to enforce the constraint
  double m_tol;                     //!< Tolerance for the constraint to be satisfied. 
  double m_rescale;                 //!< Rescale the constraint (e.g., sphere radius) by this much in total
  int m_rescale_steps;              //!< Rescale constrain over this many step
  int m_rescale_freq;               //!< Skip this many steps between rescaling constrain
  double m_scale;                   //!< Rescale the constraint (e.g., sphere radius) by this much in each step (=m_rescale**(m_rescale_freq/m_rescale_steps))
  string m_group;                   //!< Apply constraint only to particles in this group
  
};

typedef shared_ptr<Constraint> ConstraintPtr;  //!< Shared pointer to the Constraint object

#endif



