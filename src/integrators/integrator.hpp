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
 * \file integrator.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Declaration of Integrator class
 */ 

#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

#include <string>
#include <memory>


#include "system.hpp"
#include "messenger.hpp"
#include "potential.hpp"
#include "constrainer.hpp"
#include "aligner.hpp" 
#include "value.hpp"
#include "parse_parameters.hpp"

using std::string;

using std::make_shared;

/*! Integrator class is the base class for handling different numerical 
 *  integrators for integrating equations of motion (e.g., NVE).
 *  This is an abstract class and its children will implement 
 *  actual integrators.
*/
class Integrator
{
public:
  
  //! Constructor
  //! \param sys Pointer to a System object containing all particles
  //! \param msg Internal message handler
  //! \param pot Pairwise and external interaction handler
  //! \param align Pairwise and external alignment handler
  //! \param nlist Neighbour list object
  //! \param cons Enforces constraints 
  //! \param temp Temperature control object
  //! \param param Contains information about all parameters 
  Integrator(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist, ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : m_system(sys),
                                                                                                                                                                   m_msg(msg),
                                                                                                                                                                   m_potential(pot),
                                                                                                                                                                   m_align(align),
                                                                                                                                                                   m_nlist(nlist),
                                                                                                                                                                   m_constrainer(cons),
                                                                                                                                                                   m_temp(temp)
  { 
    m_known_params.push_back("dt");
    m_known_params.push_back("group");
    m_known_params.push_back("temperature_control");
    m_known_params.push_back("min_val");
    m_known_params.push_back("max_val");
    m_known_params.push_back("steps");
    if (param.find("dt") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"Integrator is setting its own time step size. While this is allowed, it may cause odd behaviour in case different integrators set different internal step size.");
      m_msg->msg(Messenger::INFO,"Internal integrator time step set to "+param["dt"]+".");
      m_dt = lexical_cast<double>(param["dt"]);
    }
    else
    {
      m_dt = m_system->get_integrator_step();
      m_msg->msg(Messenger::INFO,"Using global (system-wide) integrator step size set to "+lexical_cast<string>(m_dt)+".");
    }
    if (m_dt <= 0.0)
    {
      m_msg->msg(Messenger::ERROR,"Integrator step size has to be larger than 0. Did you forget to use \"timestep\" command?");
      throw runtime_error("Integrator step has not been set.");
    }
    m_msg->write_config("integrator.dt",lexical_cast<string>(m_dt));
    if (param.find("group") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No group has been set. Assuming group 'all'.");
      m_group_name = "all";
    }
    else
    {
      m_group_name = param["group"];
      m_msg->msg(Messenger::INFO,"Applying integrator to group "+m_group_name+".");
    }
    m_msg->write_config("integrator.group",lexical_cast<string>(m_group_name));
  }

  virtual ~Integrator() { }
  
  //! Propagate system for a time step
  virtual void integrate() = 0;
  
  //! Check if there are no illegal parameters
  string params_ok(pairs_type& params)
  {
    for (pairs_type::iterator it_p = params.begin(); it_p != params.end(); it_p++)
      if (find(m_known_params.begin(),m_known_params.end(),(*it_p).first) == m_known_params.end())
        return (*it_p).first;
    return "";
  }

protected:
  
  SystemPtr m_system;            //!< Pointer to the System object
  MessengerPtr m_msg;            //!< Pointer to the messenger object
  PotentialPtr m_potential;      //!< Pointer to the interaction handler 
  AlignerPtr m_align;            //!< Pointer to alignment handler
  NeighbourListPtr m_nlist;      //!< Pointer to the neighbour list object
  ConstrainerPtr m_constrainer;  //!< Pointer to the handler for constraints
  ValuePtr m_temp;               //!< Pointer to the handler of current value of temperature 
  double m_dt;                   //!< time step
  string m_group_name;           //!< Name of the group to apply this integrator to
  vector<string> m_known_params; //!< Lists all known parameters accepted by a given integrator

};

typedef shared_ptr<Integrator> IntegratorPtr;

#endif
