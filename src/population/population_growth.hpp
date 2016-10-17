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
 * \file population_growth.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Jun-2015
 * \brief Declaration of PopulationGrowth class.
 */ 

#ifndef __POPULATION_GROWTH_HPP__
#define __POPULATION_GROWTH_HPP__

#include <list>
#include <stdexcept>

#include "population.hpp"
#include "rng.hpp"


using std::list;
using std::runtime_error;

/*! PopulationGrowth class handles very simple particle growth (change of the particle radius)
 *
*/
class PopulationGrowth : public Population
{
public:
  
  //! Construct PopulationGrowth object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  //! \param param Reference to the parameters that control given population control
  PopulationGrowth(SystemPtr sys, const MessengerPtr msg, pairs_type& param) : Population(sys,msg,param)
  { 
    if (param.find("scale") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Growth population control. No scale for the growth set. Assuming 2.0.");
      m_rescale = 2.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Growth population control. Scale for growth set to "+param["scale"]+".");
      m_rescale = lexical_cast<double>(param["scale"]);
    }
    m_msg->write_config("population.growth.scale",lexical_cast<string>(m_rescale));
    if (param.find("steps") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Growth population control. No number of steps for the growth set. Assuming 1000.");
      m_rescale_steps = 1000;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Growth population control. Number of steps for growth set to "+param["steps"]+".");
      m_rescale_steps = lexical_cast<double>(param["steps"]);
    }
    m_msg->write_config("population.growth.steps",lexical_cast<string>(m_rescale_steps));
    m_scale = pow(m_rescale,static_cast<double>(m_freq)/static_cast<double>(m_rescale_steps));
  }
  
  //! Particle division (emulates cell division)
  void divide(int time) { }
  
  //! Remove particle (emulates cell death)
  void remove(int time) { }
  
  //! Add particle (has no direct biological application)
  void add(int t) { }
  
  //! Change particle radius
  void grow(int);
  
  //! Change particle length
  void elongate(int time) { }
  
  
private:
  
  double m_rescale;            //!< Rescale particle radius by total of this amount
  int m_rescale_steps;         //!< Rescale particle radius over this many steps
  double m_scale;              //!< Rescale particle radius by this much in each step (=m_rescale**(m_freq/m_rescale_steps))
   
};

typedef shared_ptr<PopulationGrowth> PopulationGrowthPtr;

#endif
