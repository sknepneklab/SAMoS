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
 * \file population_density.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \author Silke Henkes, shenkes@abdn.ac.uk
 * \date 11-Dec-2014
 * \brief Declaration of PopulationDensity class.
 */ 

#ifndef __POPULATION_DENSITY_HPP__
#define __POPULATION_DENSITY_HPP__

#include <list>
#include <stdexcept>

#include "population.hpp"
#include "rng.hpp"


using std::list;
using std::runtime_error;
using std::cerr;

/*! PopulationDensity class handles density dependant particle division and death. 
 *  We set division rate and death rate. Particles are divided based on their local density and removed
 *  based on the their age, by a random event.
 * 
 *  Particles are always divided along the direction of the orientation vector 
 *  n.
 *
*/
class PopulationDensity : public Population
{
public:
  
  //! Construct PopulationDensity object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  //! \param param Reference to the parameters that control given population control
  PopulationDensity(SystemPtr sys, const MessengerPtr msg, pairs_type& param) : Population(sys,msg,param)
  { 
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("population.density.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Density population control. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("population.density.seed",param["seed"]);
    }
    if (param.find("division_rate") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No division rate set. Using default 0.001.");
      m_div_rate = 0.001;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Density population control. Setting division rate "+param["division_rate"]+".");
      m_div_rate =  lexical_cast<double>(param["division_rate"]);
      if (m_div_rate < 0.0)
      {
        m_msg->msg(Messenger::ERROR,"Density population control. Division probability has to be positive.");
        throw runtime_error("Wrong division probability.");
      }
    }
    m_msg->write_config("population.density.division_rate",lexical_cast<string>(m_div_rate));
    if (param.find("death_rate") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No death rate set. Using default 0.001.");
      m_death_rate = 0.001;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Density population control. Setting death rate "+param["death_rate"]+".");
      m_death_rate =  lexical_cast<double>(param["death_rate"]);
    }
    m_msg->write_config("population.density.death_rate",lexical_cast<string>(m_death_rate));
    if (param.find("change_prob_1") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No particle type change probability for first child set. Using default 0.");
      m_type_change_prob_1 = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Density population control. Setting particle type change probability for first child to "+param["change_prob_1"]+".");
      m_type_change_prob_1 =  lexical_cast<double>(param["change_prob_1"]);
    }
    m_msg->write_config("population.density.change_prob_1",lexical_cast<string>(m_type_change_prob_1));
    if (param.find("change_prob_2") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No particle type change probability for second child set. Using default 0.");
      m_type_change_prob_2 = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Density population control. Setting particle type change probability for second child to "+param["change_prob_2"]+".");
      m_type_change_prob_2 =  lexical_cast<double>(param["change_prob_2"]);
    }    
    m_msg->write_config("population.density.change_prob_2",lexical_cast<string>(m_type_change_prob_2));
    if (param.find("new_type") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No new particle type set. Using default 0.");
      m_new_type = 0;  // No type change
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Density population control. Type of new particle type set to "+param["new_type"]+".");
      m_new_type =  lexical_cast<int>(param["new_type"]);
    }
    m_msg->write_config("population.density.new_type",lexical_cast<string>(m_new_type));
    if (param.find("new_r") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No new particle radius set. Using default 0 (i.e. inherit from mother).");
      m_new_radius = 0.0;  // No radius change: inherit from mother
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Density population control. Radius of new particle set to "+param["new_r"]+".");
      m_new_radius =  lexical_cast<double>(param["new_r"]);
    }
    m_msg->write_config("population.density.new_radius",lexical_cast<string>(m_new_radius));
    if (param.find("poly")== param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No polydispersity set. Using default 0.");
      m_poly = 0.0;  // no polydispersity (added)
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Density population control. Polydispersity of new particles set to "+param["poly"]+".");
      m_poly =  lexical_cast<double>(param["poly"]);
    }
    if (param.find("old_group") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No old group set. Using default \"all\".");
      m_old_group = "all";
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Density population control. Old group set to "+param["old_group"]+".");
      m_old_group = param["old_group"];
    }
    m_msg->write_config("population.density.old_group",m_old_group);
    if (param.find("new_group") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No new group set. Using default \"all\".");
      m_new_group = "all";
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Density population control. New group set to "+param["new_group"]+".");
      m_new_group = param["new_group"];
    }
    m_msg->write_config("population.density.new_group",m_new_group);
    if (param.find("split_distance") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No split distance set. Using default 0.5.");
      m_split_distance = 0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Density population control. Split distance set to "+param["split_distance"]+".");
      m_split_distance = lexical_cast<double>(param["split_distance"]);
    }
    m_msg->write_config("population.density.split_distance",lexical_cast<string>(m_split_distance));
    if (param.find("rho_max") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No maximum density set. Using default 8.");
      m_rho_max = 8.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Density population control. Maximum density set to "+param["rho_max"]+".");
      m_rho_max = lexical_cast<double>(param["rho_max"]);
    }
    m_msg->write_config("population.density.rho_max",lexical_cast<string>(m_rho_max));
    if (param.find("move_split") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Density population control. No split of the offspring separation set. Assuming 0.5.");
      m_alpha = 0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Split of the offspring separation set to "+param["move_split"]+".");
      m_alpha = lexical_cast<double>(param["move_split"]);
    }
    m_msg->write_config("population.density.move_split",lexical_cast<string>(m_alpha));
    
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
  
  RNGPtr  m_rng;                 //!< Density number generator
  double  m_div_rate;            //!< Rate of division
  double  m_death_rate;          //!< Rate of death
  double  m_type_change_prob_1;  //!< Probability with which the particle type can change (particle 1, that is original particle)
  double  m_type_change_prob_2;  //!< Probability with which the particle type can change (particle 2, that is new particle)
  int m_new_type;                //!< What is type of the new particle (0 no type change)
  double m_new_radius;           //!< What is radius of new particle (0.0 no radius change)
  string m_old_group;            //!< What is old group of the particle (where to change from)
  string m_new_group;            //!< What is group of new particle (what to change to)
  double m_split_distance;       //!< Fraction of the particle radius to split after the division
  double m_rho_max;              //!< Maximum local density where division rate decays 
  double m_alpha;                //!< When dividing particles, move new one to alpha*m_split_distance and the old one to (1-alpha)*split_distance
  double m_poly;				         //!< When dividing particles, give the daughter a radius of m_new_radius with polydispersity poly. To avoid selecting for small radii.
   
};

typedef shared_ptr<PopulationDensity> PopulationDensityPtr;

#endif
