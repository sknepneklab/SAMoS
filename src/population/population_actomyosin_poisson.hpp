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
 * \file population_actomyosin_poisson.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 26-Jul-2016
 * \brief Declaration of PopulationActomyosinPoisson class.
 */ 

#ifndef __POPULATION_ACTOMYOSIN_POISSON_HPP__
#define __POPULATION_ACTOMYOSIN_POISSON_HPP__

#include <list>
#include <stdexcept>
#include <cmath>

#include "population.hpp"
#include "rng.hpp"


using std::list;
using std::runtime_error;
using std::sqrt;
using std::exp;

/*! PopulationActomyosinPoisson class handles a basic model for actomyosin binding and unbiding.
 *  This model uses Poisson process based on the particle age to decide if it binds or unbinds.
 *
*/
class PopulationActomyosinPoisson : public Population
{
public:
  
  //! Construct PopulationActomyosinPoisson object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  //! \param param Reference to the parameters that control given population control
  PopulationActomyosinPoisson(SystemPtr sys, const MessengerPtr msg, pairs_type& param) : Population(sys,msg,param)
  { 
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyosin Poisson population control. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("population.actomyosin_poisson.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyosin Poisson population control. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("population.actomyosin_poisson.seed",param["seed"]);
    }
    
    if (param.find("detachment_rate") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyosin Poisson population control. No detachment rate set. Using default 0.5.");
      m_detach_rate = 0.5/(static_cast<double>(m_freq)*m_system->get_integrator_step());
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyosin Poisson population control. Setting detachment rate "+param["detachment_rate"]+".");
      m_detach_rate = lexical_cast<double>(param["detachment_rate"])/(static_cast<double>(m_freq)*m_system->get_integrator_step());  
      if (m_detach_rate < 0.0)
      {
        m_msg->msg(Messenger::ERROR,"Actomyosin Poisson population control. Detachment rate has to be greater than 0.");
        throw runtime_error("Wrong detachment rate.");
      }
    }
    m_msg->write_config("population.actomyosin_poisson.detachment_rate",lexical_cast<string>(m_detach_rate));
    
    if (param.find("attachment_rate") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyosin Poisson population control. No attachment rate set. Using default 0.5.");
      m_attach_rate = 0.5/(static_cast<double>(m_freq)*m_system->get_integrator_step());
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyosin Poisson population control. Setting attachment rate "+param["attachment_rate"]+".");
      m_attach_rate = lexical_cast<double>(param["attachment_rate"])/(static_cast<double>(m_freq)*m_system->get_integrator_step());
      if (m_attach_rate < 0.0)
      {
        m_msg->msg(Messenger::ERROR,"Actomyosin Poisson population control. Attachment rate has to be greater than 0.");
        throw runtime_error("Wrong attachment rate.");
      }
    }
    m_msg->write_config("population.actomyosin_poisson.attachment_rate",lexical_cast<string>(m_attach_rate));
    
    if (param.find("lambda") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyosin Poisson population control. No attachment lambda set. Using default 0.1.");
      m_lambda = 0.1;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyosin Poisson population control. Setting attachment lambda "+param["lambda"]+".");
      m_lambda = lexical_cast<double>(param["lambda"]);
    }
    m_msg->write_config("population.actomyosin_poisson.lambda",lexical_cast<string>(m_lambda));
    
    if (param.find("re") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyosin Poisson population control. No attachment re set. Using default 1.0.");
      m_re = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyosin Poisson population control. Setting attachment re "+param["re"]+".");
      m_re = lexical_cast<double>(param["re"]);
    }
    m_msg->write_config("population.actomyosin_poisson.re",lexical_cast<string>(m_re));
    
    if (param.find("detached_type") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyosin Poisson population control. No type of detached bead set. Using default 3.");
      m_type_d = 3;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyosin Poisson population control. Setting type of detached bead "+param["detached_type"]+".");
      m_type_d = lexical_cast<int>(param["detached_type"]);
    }
    m_msg->write_config("population.actomyosin_poisson.detached_type",lexical_cast<string>(m_type_d));
   
    if (param.find("attached_type") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyosin population control. No type of attachment bead set. Using default 4.");
      m_type_a = 4;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyosin Poisson population control. Setting type of attachment bead "+param["attached_type"]+".");
      m_type_a = lexical_cast<int>(param["attached_type"]);
    }
    m_msg->write_config("population.actomyosin_poisson.attached_type",lexical_cast<string>(m_type_a));

    if (param.find("actin_type") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyosin Poisson population control. No type of actin bead set. Using default 1.");
      m_type_actin = 1;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyosin Poisson population control. Setting type of actin bead "+param["actin_type"]+".");
      m_type_actin = lexical_cast<int>(param["actin_type"]);
    }
    m_msg->write_config("population.actomyosin_poisson.actin_type",lexical_cast<string>(m_type_actin));   
  }
  
  //! This funciton controls attachment (note: function name derives from the intial intent of the population classes to use to treat cell division)
  void divide(int);
  
  //! This function controls detachement
  void remove(int) { }
  
  //! Not used here
  void add(int t) { }
  
  //! Not used here 
  void grow(int time) { }
  
  //! Not used here
  void elongate(int time) { }
  
  
private:
  
  RNGPtr m_rng;                  //!< Actomyosin number generator
  double m_detach_rate;          //!< Bare detachment rate
  double m_lambda;               //!< Lambda parameter for detachment in the presence of force
  double m_re;                   //!< Range over which the detachement is affected by the force
  double m_attach_rate;          //!< Bare attachement rate
  int m_type_d;                  //!< Type when detached
  int m_type_a;                  //!< Type when attached
  int m_type_actin;              //!< Type of the actin beads

};

typedef shared_ptr<PopulationActomyosinPoisson> PopulationActomyosinPoissonPtr;

#endif
