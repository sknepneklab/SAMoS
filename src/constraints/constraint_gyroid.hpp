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
 *  given by equation \f$ \cos \left(\frac{2 \pi  x}{lx}\right) \sin \left(\frac{2 \pi y}{ly}\right)+\sin \left(\frac{2 \pi  x}{lx}\right) \cos
   \left(\frac{2 \pi  z}{lz}\right)+\cos \left(\frac{2 \pi y}{ly}\right) \sin \left(\frac{2 \pi  z}{lz}\right) = 0. \f$
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
    m_msg->write_config("constraint.gyroid.type","gyroid");
  }
  
  //! Computes normal to the surface
  void compute_normal(Particle&, double&, double&, double&); 
   
  // Computer gradient at a point
  void compute_gradient(Particle&, double&, double&, double&);
  
  // Value of the constraint
  double constraint_value(Particle&); 
    
  
};

typedef shared_ptr<ConstraintGyroid> ConstraintGyroidPtr;  //!< Shared pointer to the Constraint object

#endif