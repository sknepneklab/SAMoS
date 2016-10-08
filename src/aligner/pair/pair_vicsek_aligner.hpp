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
 * \file pair_vicsek_aligner.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Jan-2014
 * \brief Declaration of PairVicsekAlign class
 */ 

#ifndef __PAIR_VICSEK_ALIGNER_HPP__
#define __PAIR_VICSEK_ALIGNER_HPP__

#include <cmath>

#include "pair_aligner.hpp"

using std::make_pair;
using std::sqrt;

/*! PairVicsekAlign is an alignment rule that computes the average direction of the 
 *  velocity averaged over a certain neighbourhood of the particle.
 *  \f$ \langle \vec v_i\rangle = \sum_{\mathrm{j n.n i}} \vec v_j \f$ 
 */
class PairVicsekAlign : public PairAlign
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param param Contains information about all parameters (k)
  PairVicsekAlign(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, pairs_type& param) : PairAlign(sys, msg, nlist, param)
  {
    if (param.find("rcut") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No cutoff distance (rcut) specified for Vicsek alignment. Setting it to the global neighbour list cutoff.");
      m_rcut = m_nlist->get_cutoff();
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global cutoff distance (rcut) for Vicsek alignment is set to "+param["rcut"]+".");
      m_rcut = lexical_cast<double>(param["rcut"]);
    }
  }
                                                                                                                
  //! Set pair parameters data for pairwise interactions    
  //! \note It makes no sense to distinguish different pair parameters!
  void set_pair_parameters(pairs_type& pair_param)
  {
    m_has_pair_params = false;
  }
  
  //! Returns true since Vicsek potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes "torques"
  void compute();
  
  
private:
       
  double m_rcut;     //!<  Cutoff distance (has to be less than neighbour list cutoff)
     
};

typedef shared_ptr<PairVicsekAlign> PairVicsekAlignPtr;

#endif
