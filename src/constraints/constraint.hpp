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
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

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
  Constraint(SystemPtr sys, MessengerPtr msg, pairs_type& param) : m_system(sys), m_msg(msg), m_rescale(1.0), m_rescale_steps(1000), m_rescale_freq(10) 
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
    if (param.find("rescale") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Constraint. Rescaling constraint by factor "+param["rescale"]+".");
      m_rescale = lexical_cast<double>(param["rescale"]);
    }
    if (param.find("rescale_steps") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Constraint. Rescaling of the constraint done over "+param["rescale_steps"]+" time steps.");
      m_rescale_steps = lexical_cast<double>(param["rescale_steps"]);
    }
    if (param.find("rescale_freq") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Constraint. Rescaling constraint every "+param["rescale_freq"]+" time steps.");
      m_rescale_freq = lexical_cast<double>(param["rescale_freq"]);
    }
    m_scale = pow(m_rescale,static_cast<double>(m_rescale_freq)/static_cast<double>(m_rescale_steps));
    if (m_scale != 1.0)
      m_msg->msg(Messenger::INFO,"Constraint. Rescaling constraint with scale factor"+lexical_cast<string>(m_scale)+" per time steps.");
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
  
protected:
  
  SystemPtr  m_system;              //!< Pointer to the system object
  MessengerPtr m_msg;               //!< Handles internal messages
  int m_max_iter;                   //!< Maximum number of iterations to enforce the constraint
  double m_tol;                     //!< Tolerance for the constraint to be satisfied. 
  double m_rescale;                 //!< Rescale the constraint (e.g., sphere radius) by this much in total
  int m_rescale_steps;              //!< Rescale constrain over this many step
  int m_rescale_freq;               //!< Skip this many steps between rescaling constrain
  double m_scale;                   //!< Rescale the constraint (e.g., sphere radius) by this much in each step (=m_rescale**(m_rescale_freq/m_rescale_steps))
  
};

typedef shared_ptr<Constraint> ConstraintPtr;  //!< Shared pointer to the Constraint object

#endif



