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
    if (m_a >= m_c)
    {
      m_msg->msg(Messenger::ERROR,"Toroidal constraint. Tube radius has to be smaller than torus radius.");
      throw runtime_error("Toroidal constraint. Incompatible parameters.");
    }
    if (param.find("maxiter") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Toroidal constraint. Maximum number of iterations has not been set. Assuming 100");
      m_max_iter = 100;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Toroidal constraint. Maximum number of iterations set to "+param["maxiter"]+".");
      m_max_iter = lexical_cast<int>(param["maxiter"]);
    }
    if (param.find("tol") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Toroidal constraint. Tolerance has not been set. Assuming 1e-6.");
      m_tol = 1e-6;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Toroidal constraint. Tolerance set to "+param["tol"]+".");
      m_tol = lexical_cast<int>(param["tol"]);
    }
  }
  
  //! Computes normal to the surface
  void compute_normal(Particle&, double&, double&, double&); 
   
  //! Enforce constraint
  void enforce(Particle&);
  
  //! Rotate director around normal vector to the torus
  void rotate_director(Particle&, double);
  
  //! Rotate velocity around normal vector to the torus
  void rotate_velocity(Particle&, double);
  
  //! Project torque onto normal vector onto the torus and return rotation angle change
  double project_torque(Particle&);
    
private:
  
  double m_a;     //!< Torus tube radius
  double m_c;     //!< Torus radius
  int m_max_iter; //!< Maximum number of iterations to enforce the constraint
  double m_tol;   //!< Tolerance for the constraint to be satisfied. 
  
};

typedef shared_ptr<ConstraintTorus> ConstraintTorusPtr;  //!< Shared pointer to the Constraint object

#endif