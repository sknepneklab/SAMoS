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
 * \file population.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 28-Sept-2014
 * \brief Declaration of Population class.
 */ 

#ifndef __POPULATION_HPP__
#define __POPULATION_HPP__

#include <string>
#include <vector>
#include <memory>

#include "messenger.hpp"
#include "system.hpp"
#include "neighbour_list.hpp"

#include "parse_parameters.hpp"

using std::string;
using std::vector;

using std::make_shared;

/*! Population class is the abstract parent class that handles all population control 
 *  mechanism (division and death/removal of particles).
 *
*/
class Population
{
public:
  
  //! Construct Population object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  //! \param param Reference to the parameters that control given population control
  Population(SystemPtr sys, const MessengerPtr msg, pairs_type& param) : m_system(sys), m_msg(msg), m_has_nlist(false)
  { 
    if (param.find("group") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Population. No group has been set. Assuming group 'all'.");
      m_group_name = "all";
    }
    else
    {
      m_group_name = param["group"];
      m_msg->msg(Messenger::INFO,"Population. Applying population to group "+m_group_name+".");
    }
    if (param.find("freq") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Population. No population control frequency has been set. Assuming default 100.");
      m_freq = 100;
    }
    else
    {
      m_freq = lexical_cast<int>(param["freq"]);
      m_msg->msg(Messenger::INFO,"Population. Population control frequency set to "+param["freq"]+".");     
    }
    m_msg->write_config("population.freq",lexical_cast<string>(m_freq));
  }
  
  //! Set neighbour list. 
  //! \note Neighbour list is not set in the contructor because some models do not require it. 
  //! This is a bit clumsy but gives us more flexibility
  void set_nlist(NeighbourListPtr nlist) 
  {
    m_nlist = nlist;
    m_has_nlist = true;
  }
  
  string get_group() { return m_group_name; }
  
  //! Particle division (emulates cell division)
  virtual void divide(int) = 0;
  
  //! Remove particle (emulates cell death)
  virtual void remove(int) = 0;
  
  //! Add particle (has no direct biological application)
  virtual void add(int) = 0;
  
  //! Change particle radius
  virtual void grow(int) = 0;
  
  //! Change particle length
  virtual void elongate(int) = 0;
  
protected:
  
  SystemPtr m_system;            //!< Contains pointer to the System object
  MessengerPtr m_msg;            //!< Handles messages sent to output
  NeighbourListPtr m_nlist;         //!< Handles NeighbourList object
  string m_group_name;           //!< Name of the group to apply this population to
  int m_freq;                    //!< Frequency (in the number of time steps) with which to attempt population control
  double m_l_rescale;            //!< Rescale particle length by total of this amount
  int m_l_rescale_steps;         //!< Rescale particle length over this many steps
  double m_l_scale;              //!< Rescale particle length by this much in each step (=m_l_rescale**(m_grow_l_freq/m_l_rescale_steps))
  bool m_has_nlist;              //!< If true, population has neighbour list set
   
};

typedef shared_ptr<Population> PopulationPtr;

#endif
