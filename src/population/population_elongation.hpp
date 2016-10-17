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
 * \file population_elongation.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Jun-2015
 * \brief Declaration of PopulationElongation class.
 */ 

#ifndef __POPULATION_ELONGATION_HPP__
#define __POPULATION_ELONGATION_HPP__

#include <list>
#include <stdexcept>

#include "population.hpp"
#include "rng.hpp"


using std::list;
using std::runtime_error;

/*! PopulationElongation class handles very simple particle Elongation (change of the particle radius)
 *
*/
class PopulationElongation : public Population
{
public:
  
  //! Construct PopulationElongation object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  //! \param param Reference to the parameters that control given population control
  PopulationElongation(SystemPtr sys, const MessengerPtr msg, pairs_type& param) : Population(sys,msg,param)
  { 
    if (param.find("scale") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Elongation population control. No scale for the Elongation set. Assuming 2.0.");
      m_rescale = 2.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Elongation population control. Scale for Elongation set to "+param["scale"]+".");
      m_rescale = lexical_cast<double>(param["scale"]);
    }
    m_msg->write_config("population.Elongation.scale",lexical_cast<string>(m_rescale));
    if (param.find("steps") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Elongation population control. No number of steps for the Elongation set. Assuming 1000.");
      m_rescale_steps = 1000;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Elongation population control. Number of steps for Elongation set to "+param["steps"]+".");
      m_rescale_steps = lexical_cast<double>(param["steps"]);
    }
    m_msg->write_config("population.elongation.steps",lexical_cast<string>(m_rescale_steps));
    m_scale = pow(m_rescale,static_cast<double>(m_freq)/static_cast<double>(m_rescale_steps));
  }
  
  //! Particle division (emulates cell division)
  void divide(int time) { }
  
  //! Remove particle (emulates cell death)
  void remove(int time) { }
  
  //! Add particle (has no direct biological application)
  void add(int t) { }
  
  //! Change particle radius
  void grow(int t) { }
  
  //! Change particle length
  void elongate(int);
  
  
private:
  
  double m_rescale;            //!< Rescale particle radius by total of this amount
  int m_rescale_steps;         //!< Rescale particle radius over this many steps
  double m_scale;              //!< Rescale particle radius by this much in each step (=m_rescale**(m_freq/m_rescale_steps))
   
};

typedef shared_ptr<PopulationElongation> PopulationElongationPtr;

#endif
