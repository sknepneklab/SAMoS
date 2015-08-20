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
 * \file constraint_gaussian_bump.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Jun-2015
 * \brief Declaration of ConstraintGaussianBump class.
 */ 

#ifndef __CONSTRAINT_GAUSSIAN_BUMP_HPP__
#define __CONSTRAINT_GAUSSIAN_BUMP_HPP__

#include <cmath>

#include "system.hpp"
#include "parse_parameters.hpp"
#include "constraint.hpp"


using std::exp;


/*! Enforces all particles to be on the surface of a Gaussian bump centred in the origin.
 *  Form of the bump is \f$ z(x,y)=A \exp\left(-\frac{x^2}{a^2}\right)\exp\left(-\frac{y^2}{b^2}\right) \f$,
 *  where parameter \f$ A \f$ controls hight of the bump and parameters \f$ a \f$ and \f$ b \f$ control the
 *  width of the bump in the x and y direction, respectively.
*/
class ConstraintGaussianBump : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds (e.g., torus radius)
  ConstraintGaussianBump(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    if (param.find("A") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Gaussian bump constraint. No height set set. Assuming 1.");
      m_A = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Gaussian bump constraint. Heigh set to "+param["A"]+".");
      m_A = lexical_cast<double>(param["A"]);
    }
    m_msg->write_config("constraint.gaussian_bump.A",lexical_cast<string>(m_A));
    if (param.find("a") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Gaussian bump constraint. No width in x direction set. Assuming 1.0.");
      m_a = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Gaussian bump constraint. Width in x direction set to "+param["a"]+".");
      m_a = lexical_cast<double>(param["a"]);
      m_a *= m_a;
    }
    m_msg->write_config("constraint.gaussian_bump.a",lexical_cast<string>(m_a));
    if (param.find("b") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Gaussian bump constraint. No width in y direction set. Assuming 1.0.");
      m_b = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Gaussian bump constraint. Width in y direction set to "+param["b"]+".");
      m_b = lexical_cast<double>(param["b"]);
      m_b *= m_b;
    }
    m_msg->write_config("constraint.gaussian_bump.b",lexical_cast<string>(m_b));
  }
  
  //! Computes normal to the surface
  void compute_normal(Particle&, double&, double&, double&); 
   
  // Computer gradient at a point
  void compute_gradient(Particle&, double&, double&, double&);
  
  // Value of the constraint
  double constraint_value(Particle&); 
    
    
private:
  
  double m_A;     //!< GaussianBump height
  double m_a;     //!< GaussianBump width in x direction
  double m_b;     //!< GaussianBump width in y direction
    
};

typedef shared_ptr<ConstraintGaussianBump> ConstraintGaussianBumpPtr;  //!< Shared pointer to the Constraint object

#endif