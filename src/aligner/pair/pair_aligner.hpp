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
 * \file pair_aligner.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Jan-2014
 * \brief Declaration of PairAlign class
 */ 

#ifndef __PAIR_ALIGN_HPP__
#define __PAIR_ALIGN_HPP__

#include <map>
#include <vector>
#include <string>

#include "system.hpp"
#include "neighbour_list.hpp"

#include "parse_parameters.hpp"

using std::map;
using std::pair;
using std::string;
using std::vector;

typedef map<pair<int,int>, map<string,double> > PairAlignData;

/*! PairAlign is the base class (abstract) that handles
 *  all calls to the pair potential evaluations. Children 
 *  of this class will implement actual force fields.
 */
class PairAlign
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param param Contains information about all parameters 
  PairAlign(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, pairs_type& param) : m_system(sys), 
                                                                                          m_msg(msg),
                                                                                          m_nlist(nlist),
                                                                                          m_has_pair_params(false)
                                                                                          { }
   
  //! Get the total potential energy
  double get_potential_energy() { return m_potential_energy; } //!< \return value of the total potential energy 
   
  //! Set pair parameters data for pairwise alignment    
  virtual void set_pair_parameters(pairs_type&) = 0;
  
  //! Returns true if the specific pair alignment needs neighbour list
  virtual bool need_nlist() = 0;
  
  //! Computes alignment torques for all particles
  virtual void compute() = 0;
  
protected:
       
  SystemPtr m_system;              //!< Pointer to the System object
  MessengerPtr m_msg;              //!< Handles messages sent to output
  NeighbourListPtr m_nlist;        //!< Handles NeighbourList object
  PairAlignData m_pair_params;     //!< Handles specific parameters for a given pair
  bool m_has_pair_params;          //!< Flag that controls if pair parameters are set
  double m_potential_energy;       //!< Total potential energy (\note Need some thinking on how to properly define it)
    
};


typedef shared_ptr<PairAlign> PairAlignPtr;


#endif
