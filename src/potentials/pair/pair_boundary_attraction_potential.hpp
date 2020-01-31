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
 * \file pair_boundary_attraction_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Jun-2016
 * \brief Declaration of PairBoundaryAttractionPotential class
 */ 

#ifndef __PAIR_BOUNDARY_ATTRACTION_POTENTIAL_HPP__
#define __PAIR_BOUNDARY_ATTRACTION_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;
using std::cos;
using std::sin;

//! Structure that handles parameters for the boundary attraction pair potential
struct BoundaryAttractionParameters
{
  double epsilon;
  double rc;
  double wc;
};

/*! PairBoundaryAttractionPotential implements the attraction between boundary cells and the bulk in the active vertex tissue model.
 *  Force is only existed on the boundary cells and is modelled based on the corase-grained model of Cooke, at al. (Phys. Rev. E 72, 011506 (2005)).
 *  If the distance between cell centres is between \f$ r_c $\f and \f$ r_c + w_c $\f force is equal to 
 *  \f$ -\frac{\pi\varepsilon}{2w_c}\sin\left(\frac{\pi(r_{ij}-r_c)}{w_c}\right){\hat\vec r}_{ij} \f$,
 *  where \f$ \varepsilon $\f measures interaction strength. 
 */
class PairBoundaryAttractionPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (k)
  PairBoundaryAttractionPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    if (param.find("epsilon") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No attraction strength (epsilon) specified for boundary attraction pair potential. Setting it to 1.");
      m_epsilon = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global attraction strength (epsilon) for boundary attraction pair potential is set to "+param["epsilon"]+".");
      m_epsilon = lexical_cast<double>(param["epsilon"]);
    }
    m_msg->write_config("potential.pair.boundary_attraction.epsilon",lexical_cast<string>(m_epsilon));
    
    if (param.find("rc") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No lower cutoff (rc) specified for boundary attraction pair potential. Setting it to 1.");
      m_rc = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global lower cutoff (rc) for boundary attraction pair potential is set to "+param["rc"]+".");
      m_rc = lexical_cast<double>(param["rc"]);
    }
    m_msg->write_config("potential.pair.boundary_attraction.rc",lexical_cast<string>(m_rc));
    
    if (param.find("wc") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No width (wc) specified for boundary attraction pair potential. Setting it to 1.");
      m_wc = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global width (wc) for boundary attraction pair potential is set to "+param["wc"]+".");
      m_wc = lexical_cast<double>(param["wc"]);
    }
    m_msg->write_config("potential.pair.boundary_attraction.wc",lexical_cast<string>(m_wc));
    
    if (param.find("bulk") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Boundary attraction pair potential is computed only with cells in the bulk (no boundary-boundary interactions).");
      m_exclude_boundary = true;
      m_msg->write_config("potential.pair.boundary_attraction.bulk","true");
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Boundary attraction pair potential is computed between all neighbouring cells (including boundary-boundary interactions).");
      m_exclude_boundary = false;
      m_msg->write_config("potential.pair.boundary_attraction.bulk","false");
    }
    
    m_pair_params = new BoundaryAttractionParameters*[m_ntypes];
    for (int i = 0; i < m_ntypes; i++)
    {
      m_pair_params[i] = new BoundaryAttractionParameters[m_ntypes];
      for (int j = 0; j < m_ntypes; j++)
      {
        m_pair_params[i][j].epsilon = m_epsilon;
        m_pair_params[i][j].rc = m_rc;
        m_pair_params[i][j].wc = m_wc;
      }
    }
    
  }
  
  virtual ~PairBoundaryAttractionPotential()
  {
    for (int i = 0; i < m_ntypes; i++)
      delete [] m_pair_params[i];
    delete [] m_pair_params;
  }
                                                                                                                
  //! Set pair parameters data for pairwise interactions    
  void set_pair_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    int type_1, type_2;
    
    if (pair_param.find("type_1") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in boundary attraction potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in boundary attraction potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    if (pair_param.find("epsilon") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Boundary attraction pair potential. Setting attraction strength to "+pair_param["epsilon"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["epsilon"] = lexical_cast<double>(pair_param["epsilon"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Boundary attraction pair potential. Using default attraction strength ("+lexical_cast<string>(m_epsilon)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["epsilon"] = m_epsilon;
    }
    m_msg->write_config("potential.pair.boundary_attraction.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".push",lexical_cast<string>(param["epsilon"]));
    
    if (pair_param.find("rc") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Boundary attraction pair potential. Setting lower cutoff to "+pair_param["rc"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rc"] = lexical_cast<double>(pair_param["rc"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Boundary attraction pair potential. Using default lower cutoff ("+lexical_cast<string>(m_rc)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rc"] = m_rc;
    }
    m_msg->write_config("potential.pair.boundary_attraction.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".push",lexical_cast<string>(param["rc"]));
        
    if (pair_param.find("wc") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Boundary attraction pair potential. Setting width to "+pair_param["wc"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["wc"] = lexical_cast<double>(pair_param["wc"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Boundary attraction pair potential. Using default width ("+lexical_cast<string>(m_wc)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["wc"] = m_wc;
    }
    m_msg->write_config("potential.pair.boundary_attraction.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".push",lexical_cast<string>(param["wc"]));    
    
    m_pair_params[type_1-1][type_2-1].epsilon = param["epsilon"];
    m_pair_params[type_1-1][type_2-1].rc = param["rc"];
    m_pair_params[type_1-1][type_2-1].wc = param["wc"];
    if (type_1 != type_2)
    {
      m_pair_params[type_2-1][type_1-1].epsilon = param["epsilon"];
      m_pair_params[type_2-1][type_1-1].rc = param["rc"];
      m_pair_params[type_2-1][type_1-1].wc = param["wc"];
    }
    
    m_has_pair_params = true;
  }
  
  //! Returns true since vertex-particle potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
       
  double m_epsilon;                             //!< strength of the potential
  double m_rc;                                  //!< lower cutoff
  double m_wc;                                  //!< width of the potential  
  bool m_exclude_boundary;                      //!< If true, compute force only with bulk cells
  BoundaryAttractionParameters** m_pair_params;       //!< type specific pair parameters 
     
};

typedef shared_ptr<PairBoundaryAttractionPotential> PairBoundaryAttractionPotentialPtr;

#endif
