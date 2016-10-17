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
 * \file population_cell.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Dec-2015
 * \brief Declaration of PopulationCell class.
 */ 

#ifndef __POPULATION_CELL_HPP__
#define __POPULATION_CELL_HPP__

#include <list>
#include <stdexcept>

#include "population.hpp"
#include "rng.hpp"


using std::list;
using std::runtime_error;

/*! PopulationCell class handles a basic model for cell growth, division and death. 
 *  We set the division and death rates as well as rate of cell growth, defined as a rate at which
 *  cell's native area grows. Cells are divided based on their current area and removed
 *  based on the their age, by a random event.
 * 
 *  Particles are always divided along the direction of the orientation vector 
 *  n.
 *
*/
class PopulationCell : public Population
{
public:
  
  //! Construct PopulationCell object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  //! \param param Reference to the parameters that control given population control
  PopulationCell(SystemPtr sys, const MessengerPtr msg, pairs_type& param) : Population(sys,msg,param)
  { 
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Cell population control. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("population.cell.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cell population control. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("population.cell.seed",param["seed"]);
    }
    if (param.find("division_rate") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Cell population control. No division multiplier rate set. Using default 1.0.");
      m_div_rate = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cell population control. Setting division rate "+param["division_rate"]+".");
      m_div_rate =  lexical_cast<double>(param["division_rate"]);
      if (m_div_rate <= 0.0)
      {
        m_msg->msg(Messenger::ERROR,"Cell population control. Division probability has to be greater than 0.");
        throw runtime_error("Wrong division rate.");
      }
    }
    m_msg->write_config("population.cell.division_rate",lexical_cast<string>(m_div_rate));
    if (param.find("death_rate") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Cell population control. No death rate set. Using default 0.001.");
      m_death_rate = 0.001;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cell population control. Setting death rate "+param["death_rate"]+".");
      m_death_rate =  lexical_cast<double>(param["death_rate"]);
    }
    m_msg->write_config("population.cell.death_rate",lexical_cast<string>(m_death_rate));
    if (param.find("split_distance") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Cell population control. No split distance set. Using default 0.5.");
      m_split_distance = 0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cell population control. Split distance set to "+param["split_distance"]+".");
      m_split_distance = lexical_cast<double>(param["split_distance"]);
    }
    m_msg->write_config("population.cell.split_distance",lexical_cast<string>(m_split_distance));
    if (param.find("max_area") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Cell population control. No maximum cell area set. Using default 6.");
      m_max_A0 = 6.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cell population control. Maximum native cell area set to "+param["max_area"]+".");
      m_max_A0 = lexical_cast<double>(param["max_area"]);
    }
    m_msg->write_config("population.cell.max_A0",lexical_cast<string>(m_max_A0));
    if (param.find("max_age") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Cell population control. No maximum cell age set. Using default 10,000.");
      m_max_age = 1e4;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cell population control. Maximum cell age set to "+param["max_age"]+".");
      m_max_age = lexical_cast<double>(param["max_age"]);
    }
    m_msg->write_config("population.cell.max_age",lexical_cast<string>(m_max_age));
    if (param.find("move_split") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Cell population control. No split of the offspring separation set. Assuming 0.5.");
      m_alpha = 0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cell population control. Split of the offspring separation set to "+param["move_split"]+".");
      m_alpha = lexical_cast<double>(param["move_split"]);
    }
    m_msg->write_config("population.cell.move_split",lexical_cast<string>(m_alpha));
    if (param.find("growth_rate") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Cell population control. No growth rate set. Assuming 0.01.");
      m_growth_rate = 1e-2;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cell population control. Growth rate set to "+param["growth_rate"]+".");
      m_growth_rate = lexical_cast<double>(param["growth_rate"]);
    }
    m_msg->write_config("population.cell.growth_rate",lexical_cast<string>(m_growth_rate));
    if (param.find("growth_prob") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Cell population control. No growth proability set. Assuming 0.5 (so as to introduce stochasticity)");
      m_growth_prob = 0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cell population control. Growth proability set to "+param["growth_prob"]+".");
      m_growth_prob = lexical_cast<double>(param["growth_prob"]);
    }
    m_msg->write_config("population.cell.growth_prob",lexical_cast<string>(m_growth_prob));
    if (param.find("rescale_contacts") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Cell population control. Contact distance and neigbour list distance will not be rescaled.");
      m_rescale_contacts = false;
      m_msg->write_config("population.cell.rescale_contacts","false");
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cell population control. Contact distance and neigbour list distance will be rescaled in each growth step.");
      m_rescale_contacts = true;
      m_msg->write_config("population.cell.rescale_contacts","true");
    }
    
  }
  
  //! Particle division (emulates cell division)
  void divide(int);
  
  //! Remove particle (emulates cell death)
  void remove(int);
  
  //! Add particle (has no direct biological application)
  void add(int t) { }
  
  //! Change cell native area
  void grow(int time);
  
  //! Change particle length ( makes no sense here)
  void elongate(int time) { }
  
  
private:
  
  RNGPtr m_rng;                  //!< Cell number generator
  double m_div_rate;             //!< Rate of division
  double m_death_rate;           //!< Rate of death
  double m_split_distance;       //!< Fraction of the particle radius to split after the division
  double m_max_A0;               //!< Maximum value of the cell area 
  double m_max_age;              //!< Maximum cell age
  double m_alpha;                //!< When dividing particles, move new one to alpha*m_split_distance and the old one to (1-alpha)*split_distance
  double m_growth_rate;          //!< growth (scaling) factor for the cell native area in each time step
  double m_growth_prob;          //!< probability that a cell grows in a given time step
  bool m_rescale_contacts;       //!< If true, resacle neigbour list contacnt distance as well as contact distance
  
};

typedef shared_ptr<PopulationCell> PopulationCellPtr;

#endif
