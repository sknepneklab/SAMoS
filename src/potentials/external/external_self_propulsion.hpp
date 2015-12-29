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
 * \file external_self_propulsion.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Dec-2015
 * \brief Declaration of ExternalSelfPropulsion class
 */ 

#ifndef __EXTERNAL_SELF_PROPULSION_HPP__
#define __EXTERNAL_SELF_PROPULSION_HPP__

#include "external_potential.hpp"


/*! ExternalSelfPropulsion handles active motion
 *  where each particle receives a force equal to
 *  \f$ \vec F_i = \alpha \vec n_i \f$, where \f$ \alpha \f$ 
 *  is the activity and \f$ \vec n_i \f$ is the director (polarisation). 
 *  \note This is not a potential, but to avoid code bloating we place it
 *  in the same class structure as regular external potential.
*/
class ExternalSelfPropulsion : public ExternalPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters (gravity strength)
  ExternalSelfPropulsion(SystemPtr sys, MessengerPtr msg, pairs_type& param) : ExternalPotential(sys, msg, param)  
  {
    if (param.find("alpha") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No activity specified for self propulsion. Setting it to 1.");
      m_alpha = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Activity in self propulsion set to "+param["alpha"]+".");
      m_alpha = lexical_cast<double>(param["alpha"]);
    }
    m_msg->write_config("potential.external.self_propulsion.alpha",lexical_cast<string>(m_alpha));
  }
                                                                                                                
  //! Get the total potential energy (this is not a potential so remove 0!)
  double get_potential_energy() { return 0.0; } //!< \return value of the total potential energy
  
  //! Set pair parameters data for pairwise interactions    
  void set_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    
    if (pair_param.find("type") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type has not been defined for self propulsion.");
      throw runtime_error("Missing key for self propulsion.");
    }
    
    int type = lexical_cast<int>(pair_param["type"]);
        
    if (pair_param.find("alpha") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Self propulsion. Setting alpha to "+pair_param["alpha"]+" for particle pair of type "+lexical_cast<string>(type)+".");
      param["alpha"] = lexical_cast<double>(pair_param["alpha"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Self propulsion. Using default alpha ("+lexical_cast<string>(m_alpha)+") for particle pair of type "+lexical_cast<string>(type)+").");
      param["alpha"] = m_alpha;
    }
    m_msg->write_config("potential.external.self_propulsion.type_"+pair_param["type"]+".alpha",lexical_cast<string>(param["alpha"]));

    m_type_params[type]["alpha"] = param["alpha"];
    
    m_has_params = true;
  }
  
  //! Computes forces for all particles
  void compute();
  
  
private:
       
  double m_alpha;    //!< activity
  
};

typedef shared_ptr<ExternalSelfPropulsion> ExternalSelfPropulsionPtr;

#endif
