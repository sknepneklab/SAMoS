/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
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
  Constraint(SystemPtr sys, MessengerPtr msg, pairs_type& param) : m_system(sys), m_msg(msg) 
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
  
protected:
  
  SystemPtr  m_system;              //!< Pointer to the system object
  MessengerPtr m_msg;               //!< Handles internal messages
  int m_max_iter; //!< Maximum number of iterations to enforce the constraint
  double m_tol;   //!< Tolerance for the constraint to be satisfied. 
  
};

typedef shared_ptr<Constraint> ConstraintPtr;  //!< Shared pointer to the Constraint object

#endif



