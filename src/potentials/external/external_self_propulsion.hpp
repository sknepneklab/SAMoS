/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

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

    m_type_params[type-1]["alpha"] = param["alpha"];
    
    m_has_params = true;
  }
  
  //! Computes forces for all particles
  void compute();
  
  
private:
       
  double m_alpha;    //!< activity
  
};

typedef shared_ptr<ExternalSelfPropulsion> ExternalSelfPropulsionPtr;

#endif
