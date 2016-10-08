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
 * \file population_random.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 28-Sept-2014
 * \brief Declaration of PopulationRandom class.
 */ 

#ifndef __POPULATION_RANDOM_HPP__
#define __POPULATION_RANDOM_HPP__

#include <list>
#include <stdexcept>

#include "population.hpp"
#include "rng.hpp"


using std::list;
using std::runtime_error;

/*! PopulationRandom class handles very simple particle division and death.
 *  We set division rate and death rate. Particles are divided and removed
 *  based on the their age, by a random event.
 * 
 *  Particles are always divided along the direction of the orientation vector 
 *  n.
 *
*/
class PopulationRandom : public Population
{
public:
  
  //! Construct PopulationRandom object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  //! \param param Reference to the parameters that control given population control
  PopulationRandom(SystemPtr sys, const MessengerPtr msg, pairs_type& param) : Population(sys,msg,param)
  { 
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("population.random.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("population.random.seed",param["seed"]);
    }
    if (param.find("division_rate") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No division rate set. Using default 0.001.");
      m_div_rate = 0.001;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Setting division rate "+param["division_rate"]+".");
      m_div_rate =  lexical_cast<double>(param["division_rate"]);
    }
    m_msg->write_config("population.random.division_rate",lexical_cast<string>(m_div_rate));
    if (param.find("death_rate") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No death rate set. Using default 0.001.");
      m_death_rate = 0.001;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Setting death rate "+param["death_rate"]+".");
      m_death_rate =  lexical_cast<double>(param["death_rate"]);
    }
    m_msg->write_config("population.random.death_rate",lexical_cast<string>(m_death_rate));
    if (param.find("change_prob_1") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No particle type change probability for first child set. Using default 0.");
      m_type_change_prob_1 = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Setting particle type change probability for first child to "+param["change_prob_1"]+".");
      m_type_change_prob_1 =  lexical_cast<double>(param["change_prob_1"]);
    }
    m_msg->write_config("population.random.change_prob_1",lexical_cast<string>(m_type_change_prob_1));
    if (param.find("change_prob_2") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No particle type change probability for second child set. Using default 0.");
      m_type_change_prob_2 = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Setting particle type change probability for second child to "+param["change_prob_2"]+".");
      m_type_change_prob_2 =  lexical_cast<double>(param["change_prob_2"]);
    }    
    m_msg->write_config("population.random.change_prob_2",lexical_cast<string>(m_type_change_prob_2));
    if (param.find("new_type") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No new particle type set. Using default 0.");
      m_new_type = 0;  // No type change
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Type of new particle type set to "+param["new_type"]+".");
      m_new_type =  lexical_cast<int>(param["new_type"]);
    }
    m_msg->write_config("population.random.new_type",lexical_cast<string>(m_new_type));
    if (param.find("new_r") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No new particle radius set. Using default 0.");
      m_new_radius = 0.0;  // No type change
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Radius of new particle set to "+param["new_r"]+".");
      m_new_radius =  lexical_cast<double>(param["new_r"]);
    }
    m_msg->write_config("population.random.new_r",lexical_cast<string>(m_new_radius));
    if (param.find("old_group") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No old group set. Using default \"all\".");
      m_old_group = "all";
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Old group set to "+param["old_group"]+".");
      m_old_group = param["old_group"];
    }
    m_msg->write_config("population.random.old_group",m_old_group);
    if (param.find("new_group") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No new group set. Using default \"all\".");
      m_new_group = "all";
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. New group set to "+param["new_group"]+".");
      m_new_group = param["new_group"];
    }
    m_msg->write_config("population.random.new_group",m_new_group);
    if (param.find("move_split") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No split of the offspring separation set. Assuming 0.5.");
      m_alpha = 0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Split of the offspring separation set to "+param["move_split"]+".");
      m_alpha = lexical_cast<double>(param["move_split"]);
    }
    m_msg->write_config("population.random.move_split",lexical_cast<string>(m_alpha));
    
  }
  
  //! Particle division (emulates cell division)
  void divide(int);
  
  //! Remove particle (emulates cell death)
  void remove(int);
  
  //! Add particle (has no direct biological application)
  void add(int t) { }
  
  //! Change particle radius
  void grow(int time) { }
  
  //! Change particle length
  void elongate(int time) { }
  
  
private:
  
  RNGPtr  m_rng;                 //!< Random number generator
  double  m_div_rate;            //!< Rate of division
  double  m_death_rate;          //!< Rate of death
  double  m_type_change_prob_1;  //!< Probability with which the particle type can change (particle 1, that is original particle)
  double  m_type_change_prob_2;  //!< Probability with which the particle type can change (particle 2, that is new particle)
  int m_new_type;                //!< What is type of the new particle (0 no type change)
  double m_new_radius;           //!< What is radius of new particle (0.0 no radius change)
  string m_old_group;            //!< What is old group of the particle (where to change from)
  string m_new_group;            //!< What is group of new particle (what to change to)
  double m_alpha;                //!< When dividing particles, move new one to alpha*radius and the old one to (1-alpha)*radius
   
};

typedef shared_ptr<PopulationRandom> PopulationRandomPtr;

#endif
