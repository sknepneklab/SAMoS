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
 *  given by equation \f$ \cos \left(\frac{2 \pi  x}{\text{lx}}\right) \sin \left(\frac{2 \pi y}{\text{ly}}\right)+\sin \left(\frac{2 \pi  x}{\text{lx}}\right) \cos
   \left(\frac{2 \pi  z}{\text{lz}}\right)+\cos \left(\frac{2 \pi y}{\text{ly}}\right) \sin \left(\frac{2 \pi  z}{\text{lz}}\right) = 0. \f$
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
    if (param.find("maxiter") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Gyroid constraint. Maximum number of iterations has not been set. Assuming 100");
      m_max_iter = 100;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Gyroid constraint. Maximum number of iterations set to "+param["maxiter"]+".");
      m_max_iter = lexical_cast<int>(param["maxiter"]);
    }
    if (param.find("tol") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Gyroid constraint. Tolerance has not been set. Assuming 1e-6.");
      m_tol = 1e-6;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Gyroid constraint. Tolerance set to "+param["tol"]+".");
      m_tol = lexical_cast<int>(param["tol"]);
    }
  }
  
  //! Computes normal to the surface
  void compute_normal(Particle&, double&, double&, double&); 
   
  //! Enforce constraint
  void enforce(Particle&);
  
  //! Rotate director around normal vector to the gyroid
  void rotate_director(Particle&, double);
  
  //! Rotate velocity around normal vector to the gyroid
  void rotate_velocity(Particle&, double);
  
  //! Project torque onto normal vector onto the gyroid and return rotation angle change
  double project_torque(Particle&);
    
private:
  
  int m_max_iter; //!< Maximum number of iterations to enforce the constraint
  double m_tol;   //!< Tolerance for the constraint to be satisfied. 
  
};

typedef shared_ptr<ConstraintGyroid> ConstraintGyroidPtr;  //!< Shared pointer to the Constraint object

#endif