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
class ExternalGravityPotential : public ExternalPotential
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
      m_g = lexical_cast<double>(param["g"]);
    }
    m_msg->write_config("potential.external.gravity.g",lexical_cast<string>(m_g));
  }
                                                                                                                
  //! Get the total potential energy
  double get_potential_energy() { return m_potential_energy; } //!< \return value of the total potential energy
  
  //! Set pair parameters data for pairwise interactions    
  void set_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    
    if (pair_param.find("type") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type has not been defined for external gravitational potential.");
      throw runtime_error("Missing key for external potential parameters.");
    }
    
    int type = lexical_cast<int>(pair_param["type"]);
        
    if (pair_param.find("g") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"External gravitational potential. Setting g to "+pair_param["g"]+" for particle pair of type "+lexical_cast<string>(type)+".");
      param["g"] = lexical_cast<double>(pair_param["g"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External gravitational potential. Using default g ("+lexical_cast<string>(m_g)+") for particle pair of type "+lexical_cast<string>(type)+").");
      param["g"] = m_g;
    }
    m_msg->write_config("potential.external.gravity.type_"+pair_param["type"]+".g",lexical_cast<string>(param["g"]));

    m_type_params[type-1]["g"] = param["g"];
    
    m_has_params = true;
  }
  
  //! Computes potentials and forces for all particles
  void compute();
  
  
private:
       
  double m_g;    //!< strength of external gravitational potential
  
};

typedef shared_ptr<ExternalGravityPotential> ExternalGravityPotentialPtr;

#endif
