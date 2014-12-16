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
    if (param.find("maxiter") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Ellipsoid constraint. Maximum number of iterations has not been set. Assuming 100");
      m_max_iter = 100;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Ellipsoid constraint. Maximum number of iterations set to "+param["maxiter"]+".");
      m_max_iter = lexical_cast<int>(param["maxiter"]);
    }
    if (param.find("tol") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Ellipsoid constraint. Tolerance has not been set. Assuming 1e-6.");
      m_tol = 1e-6;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Ellipsoid constraint. Tolerance set to "+param["tol"]+".");
      m_tol = lexical_cast<int>(param["tol"]);
    }
  }
  
  //! Computes normal to the surface
  void compute_normal(Particle&, double&, double&, double&); 
   
  //! Enforce constraint
  void enforce(Particle&);
  
  //! Rotate director around normal vector to the ellipsoid
  void rotate_director(Particle&, double);
  
  //! Rotate velocity around normal vector to the ellipsoid
  void rotate_velocity(Particle&, double);
  
  //! Project torque onto normal vector onto the ellipsoid and return rotation angle change
  double project_torque(Particle&);
    
private:
  
  double m_a;     //!< Ellipsoid parameter a
  double m_b;     //!< Ellipsoid parameter b
  double m_c;     //!< Ellipsoid parameter c
  int m_max_iter; //!< Maximum number of iterations to enforce the constraint
  double m_tol;   //!< Tolerance for the constraint to be satisfied. 
  
};

typedef shared_ptr<ConstraintEllipsoid> ConstraintEllipsoidPtr;  //!< Shared pointer to the Constraint object

#endif