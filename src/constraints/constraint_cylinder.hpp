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
 *   (c) 2013, 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file constraint_cylinder.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 10-Oct-2014
 * \brief Declaration of ConstraintCylinder class.
 */ 

#ifndef __CONSTRAINT_CYLINDER_HPP__
#define __CONSTRAINT_CYLINDER_HPP__

#include <cmath>

#include "system.hpp"
#include "parse_parameters.hpp"
#include "constraint.hpp"


using std::sqrt;
using std::sin;
using std::cos;

/*! Enforces all particles to be on the surface of a cylinder of radius 
 *  R along z axis passing through its centre. All velocities will point
 *  in the tangent direction.
*/
class ConstraintCylinder : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds (e.g., sphere radius)
  ConstraintCylinder(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    if (param.find("r") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Cylindrical constraint. No radius set. Assuming 1.");
      m_r = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cylindrical constraint. Radius set to "+param["r"]+".");
      m_r = lexical_cast<double>(param["r"]);
    }
  }
  
  //! Enforce constraint
  void enforce(Particle&);
  
  //! Computes normal to the surface
  void compute_normal(Particle& p, double& Nx, double& Ny, double& Nz)
  { 
    Nx = p.x/m_r; Ny = p.y/m_r; Nz = 0.0;
  }
  
  // Computer gradient at a point
  void compute_gradient(Particle& p, double& gx, double& gy, double& gz) { }
  
  // Value of the constraint
  double constraint_value(Particle& p) { return 0.0; }
  
  // Rescale constraint
  bool rescale();
      
private:
  
  double m_r;     //!< Radius of the confining cylinder
  
};

typedef shared_ptr<ConstraintCylinder> ConstraintCylinderPtr;  //!< Shared pointer to the Constraint object

#endif