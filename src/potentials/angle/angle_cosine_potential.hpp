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
 * \file angle_cosine_potential.hpp
 * \author Amit Das, ncbs, India, dosamit@gmail.com
 * \date 05-Nov-2014
 * \brief Declaration of AngleCosinePotential class
 */ 

#ifndef __ANGLE_COSINE_POTENTIAL_HPP__
#define __ANGLE_COSINE_POTENTIAL_HPP__

#include <cmath>

#include "angle_potential.hpp"

using std::make_pair;
using std::sqrt;
using std::acos;

//! Structure that handles parameters for the cosine angle
struct AngleCosineParameters
{
  double k;    // stiffness
};

/*! AngleCosinePotential implements standard cosine angle potential 
 *  Potential is given as \f$ U_{cos}\left(\theta\right) = k_\theta\left(1+cos\theta\right) \f$,
 *  where \f$ k_\theta \f$ is the stiffness, \f$ \theta \f$ is the angle between two bonds originating from the particle
 *  .
 */
class AngleCosinePotential : public AnglePotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters (k and t0)
  AngleCosinePotential(SystemPtr sys, MessengerPtr msg,  pairs_type& param) : AnglePotential(sys, msg, param)
  {
    int ntypes = m_system->get_n_angle_types();
    if (param.find("k") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No stiffness (k) specified for cosine angle potential. Setting it to 1.");
      m_k = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global stiffness (k) for cosine angle potential is set to "+param["k"]+".");
      m_k = lexical_cast<double>(param["k"]);
    }
    m_msg->write_config("potential.angle.cosine.k",lexical_cast<string>(m_k));
    m_angle_params = new AngleCosineParameters[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
      m_angle_params[i].k = m_k;
    }
  }
  
  virtual ~AngleCosinePotential()
  {
    delete [] m_angle_params;
  }
                                                                                                                
  //! Set parameters data for specific angle type    
  void set_angle_parameters(pairs_type& angle_param)
  {
    map<string,double> param;
    
    int type;
    
    if (angle_param.find("type") == angle_param.end())
    {
      m_msg->msg(Messenger::ERROR,"Angle type has not been defined for cosine angle potential.");
      throw runtime_error("Missing key for cosine angle parameters.");
    }
    type = lexical_cast<int>(angle_param["type"]);
        
    if (angle_param.find("k") != angle_param.end())
    {
      m_msg->msg(Messenger::INFO,"Cosine angle potential. Setting stiffness to "+angle_param["k"]+" for angle of types "+lexical_cast<string>(type)+".");
      param["k"] = lexical_cast<double>(angle_param["k"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cosine angle potential. Using default stiffness ("+lexical_cast<string>(m_k)+") for angles of types "+lexical_cast<string>(type)+".");
      param["k"] = m_k;
    }
    m_msg->write_config("potential.angle.cosine.type_"+angle_param["type"]+".k",lexical_cast<string>(param["k"]));
        
    m_angle_params[type-1].k = param["k"];
         
    m_has_angle_params = true;
  }
  
   
  //! Computes potentials and forces for all particles
  void compute();
  
  
private:
        
  double m_k;       //!< stiffness
  AngleCosineParameters* m_angle_params;   //!< type specific angle parameters 
    
};

typedef shared_ptr<AngleCosinePotential> AngleCosinePotentialPtr;

#endif
