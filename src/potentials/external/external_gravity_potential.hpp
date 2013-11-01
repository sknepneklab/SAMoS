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
 * \file external_gravity_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 23-Oct-2013
 * \brief Declaration of ExternalGravityPotential class
 */ 

#ifndef __EXTERNAL_GRAVITY_POTENTIAL_HPP__
#define __EXTERNAL_GRAVITY_POTENTIAL_HPP__

#include "external_potential.hpp"


/*! ExternalGravityPotential handles external gravitational
 *  potential that is for simplicity assumed to act in the z-direction
 *  (downwards). In other words,
 *  \f$ V_{grav}\left(r\right) = g z \f$ and the force is \f$ \vec F = - m g \vec e_z \f$,
 *  where \f$ g \f$ is gravitational acceleration.
*/
class ExternalGravityPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters (gravity strength)
  ExternalGravityPotential(SystemPtr sys, MessengerPtr msg, pairs_type& param) : ExternalPotential(sys, msg, param)  
  {
    if (param.find("g") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No strength specified for the external gravitational potential. Setting it to 1.");
      m_g = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Strength of the gravitational potential set to "+param["g"]+".");
      m_g = lexical_cast<double>param["g"];
    }
  }
                                                                                                                
  //! Get the total potential energy
  double get_potential_energy() { return m_potential_energy; } //!< \return value of the total potential energy
  
  //! Set pair parameters data for pairwise interactions    
  void set_parameters(int type, pairs_type& pair_param)
  {
    map<string,double> param;
    
    if (pair_param.find("g") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"External gravitational potential. Setting g to "+pair_param["g"]+" for particle pair of type "+lexical_cast<string>(type)+".");
      param["g"] = lexical_cast<double>pair_param["g"];
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External gravitational potential. Using default g ("+lexical_cast<string>(m_g)+") for particle pair of type "+lexical_cast<string>(type)+").");
      param["g"] = m_g;
    }

    m_type_params[type]["g"] = param["g"];
    
    m_has_params = true;
  }
  
  //! Computes potentials and forces for all particles
  void compute();
  
  
private:
       
  double m_g;    //!< strength of external gravitational potential
  
};

typedef shared_ptr<ExternalGravityPotential> ExternalGravityPotentialPtr;

#endif
