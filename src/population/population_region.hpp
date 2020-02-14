/* ***************************************************************************
 *
 *  Copyright (C) 2013-2020 University of Dundee
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
 * \file population_region.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 114-Feb-2020
 * \brief Declaration of PopulationRegion class.
 */ 

#ifndef __POPULATION_REGION_HPP__
#define __POPULATION_REGION_HPP__

#include <list>
#include <stdexcept>

#include "population.hpp"


using std::list;
using std::runtime_error;

/*! PopulationRegion class handles removal of cells that are outside a given region of space.
 *  If coordinates of a cell fall outside the "allowed" region, the cell is marked for removal a
 *  removed from the system. Boundaries are updated accordingly.
 * 
*/
class PopulationRegion : public Population
{
public:
  
  //! Construct PopulationRegion object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  //! \param param Reference to the parameters that control given population control
  PopulationRegion(SystemPtr sys, const MessengerPtr msg, pairs_type& param) : Population(sys,msg,param),
                                                                               m_has_left_bound{false},
                                                                               m_has_right_bound{false},
                                                                               m_has_top_bound{false},
                                                                               m_has_bottom_bound{false},
                                                                               m_xmin{0.0},
                                                                               m_xmax{0.0},
                                                                               m_ymin{0.0},
                                                                               m_ymax{0.0}
  { 
    if (param.find("xmin") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Region population control. No xmin set. Ignoring left limit.");
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Region population control. Left boundary set to "+param["xmin"]+".");
      m_xmin =  lexical_cast<double>(param["xmin"]);
      m_has_left_bound = true;
      m_msg->write_config("population.region.xmin",lexical_cast<string>(m_xmin));
    }
    if (param.find("xmax") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Region population control. No xmax set. Ignoring right limit.");
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Region population control. Right boundary set to "+param["xmax"]+".");
      m_xmax =  lexical_cast<double>(param["xmax"]);
      m_has_right_bound = true;
      m_msg->write_config("population.region.xmax",lexical_cast<string>(m_xmax));
    }
    if (param.find("ymin") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Region population control. No ymin set. Ignoring bottom limit.");
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Region population control. Bottom boundary set to "+param["ymin"]+".");
      m_ymin =  lexical_cast<double>(param["ymin"]);
      m_has_bottom_bound = true;
      m_msg->write_config("population.region.ymin",lexical_cast<string>(m_ymin));
    }
    if (param.find("ymax") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Region population control. No ymax set. Ignoring top limit.");
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Region population control. Top boundary set to "+param["ymax"]+".");
      m_ymax =  lexical_cast<double>(param["ymax"]);
      m_has_top_bound = true;
      m_msg->write_config("population.region.ymax",lexical_cast<string>(m_ymax));
    }
    
  }
  
  //! Particle division (emulates Region division)
  void divide(int time) { }
  
  //! Remove particle (emulates Region death)
  void remove(int);
  
  //! Add particle (has no direct biological application)
  void add(int t) { }
  
  //! Change Region native area
  void grow(int time) { }
  
  //! Change particle length ( makes no sense here)
  void elongate(int time)  {  }

private:

  bool m_has_left_bound;    //!< If true, there is a left boundary beyond which cells are removed
  bool m_has_right_bound;   //!< If true, there is a right boundary beyond which cells are removed
  bool m_has_top_bound;     //!< If true, there is a top boundary beyond which cells are removed
  bool m_has_bottom_bound;  //!< If true, there is a bottom boundary beyond which cells are removed
  double m_xmin;            //!< Remove all cells that have x coordinates less than this value
  double m_xmax;            //!< Remove all cells that have x coordinates greater than this value
  double m_ymin;            //!< Remove all cells that have y coordinates less than this value
  double m_ymax;            //!< Remove all cells that have y coordinates greater than this value
};

typedef shared_ptr<PopulationRegion> PopulationRegionPtr;

#endif
