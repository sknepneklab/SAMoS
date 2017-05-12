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
 * \file pair_active_nematic_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 05-May-2017
 * \brief Declaration of PairActiveNematicPotential class
 */ 

#ifndef __PAIR_ACTIVE_NEMATIC_POTENTIAL_HPP__
#define __PAIR_ACTIVE_NEMATIC_POTENTIAL_HPP__

#include "pair_potential.hpp"


//! Structure that handles parameters for the active nematic pair potential
struct ActiveNematicParameters
{
  double alpha;
  double rcut;
};

/*! PairActiveNematicPotential implements the model for active motion in active nematics
 *  introduced by F. Alaimo, C. Köhler, Axel Voigt, arXiv:1703.03707v1
 *
 *  Active force is introduced as /f$ \vec F_{active} = -\alpha \sum_j \hat{Q}_j \frac{\vec r_i - \vec r_j}{\left|\vec r_i - \vec r_j\right|^2} $\f,
 *  where \f$ \alpha \f$ is hte activity, $\hat{Q}_j^{\beta\gamma} = \left(n_j^\beta n_j^\gamma\frac{\delta^{\beta\gamma}}{3}\right) $\f 
 *  is the nematic order parameter, \f$ \vec r_i \f$ is the position of particle \f$ i \f$ and \f$ \vec n_i \f$ is its director. 
 *  Indices $\beta,\gamma\in\{x,y,z\}$.
 *
 *  \note: As many models under pair directory, this is really not a potential.
 *  \todo: Rename these directories to reflect the fact that many models are not really potentials. 
 */
class PairActiveNematicPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (alpha, sigma, and rcut)
  PairActiveNematicPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    if (!m_nlist)
    {
      m_msg->msg(Messenger::ERROR,"Active nematic pair potential requires neighbour list. None given.");
      throw runtime_error("Neighbour list required by active nematic potential, but non specified.");
    }
    if (param.find("alpha") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No contractile/extensile stress (alpha) specified for the active nematic pair potential. Setting it to 1.");
      m_alpha = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global contractile/extensile stress (alpha) for active nematic pair potential is set to "+param["alpha"]+".");
      m_alpha = lexical_cast<double>(param["alpha"]);
    }
    m_msg->write_config("potential.pair.active_nematic.alpha",lexical_cast<string>(m_alpha));
    
    if (param.find("rcut") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No cutoff distance (rcut) specified for the active nematic pair potential. Setting it to the neighbour list cutoff.");
      m_rcut = m_nlist->get_cutoff();
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global cutoff distance (rcut) for active nematic pair potential is set to "+param["rcut"]+".");
      m_rcut = lexical_cast<double>(param["rcut"]);
    }
    m_msg->write_config("potential.pair.active_nematic.rcut",lexical_cast<string>(m_rcut));


    if (m_rcut > m_nlist->get_cutoff())
      m_msg->msg(Messenger::WARNING,"Neighbour list cutoff distance (" + lexical_cast<string>(m_nlist->get_cutoff())+
      " is smaller than the active nematic cuttoff distance ("+lexical_cast<string>(m_rcut)+
      "). Results will not be reliable.");
    
    m_msg->write_config("potential.pair.active_nematic.rcut",lexical_cast<string>(m_rcut));
    
    m_pair_params = new ActiveNematicParameters*[m_ntypes];
    for (int i = 0; i < m_ntypes; i++)
    {
      m_pair_params[i] = new ActiveNematicParameters[m_ntypes];
      for (int j = 0; j < m_ntypes; j++)
      {
        m_pair_params[i][j].alpha = m_alpha;
      }
    }
  }
  
  virtual ~PairActiveNematicPotential()
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
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in active nematic potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in active nematic potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    
    if (pair_param.find("alpha") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Active nematic pair potential. Setting contractile/active stress to "+pair_param["alpha"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["alpha"] = lexical_cast<double>(pair_param["alpha"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"active nematic pair potential. Using default alpha ("+lexical_cast<string>(m_alpha)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["alpha"] = m_alpha;
    }
    if (pair_param.find("rcut") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Active nematic pair potential. Setting rcut to "+pair_param["rcut"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = lexical_cast<double>(pair_param["rcut"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Active nematic pair potential. Using default rcut ("+lexical_cast<string>(m_rcut)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = m_rcut;
    }
    m_msg->write_config("potential.pair.active_nematic.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".alpha",lexical_cast<string>(param["alpha"]));
    m_msg->write_config("potential.pair.active_nematic.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".rcut",lexical_cast<string>(param["rcut"]));   
    
    m_pair_params[type_1-1][type_2-1].alpha = param["alpha"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].alpha = param["alpha"];
    m_pair_params[type_1-1][type_2-1].rcut = param["rcut"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].rcut = param["rcut"];
    
    if (param["rcut"] > m_nlist->get_cutoff())
      m_msg->msg(Messenger::WARNING,"Neighbour list cutoff distance (" + lexical_cast<string>(m_nlist->get_cutoff())+
      "is smaller than the active nematic cuttof distance ("+lexical_cast<string>(m_rcut)+
      ") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+"). Results will not be reliable.");
    
    m_has_pair_params = true;
  }
  
  //! Returns true since active potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
       
  double m_alpha;  //!< contractile/extensile stress strength 
  double m_rcut;   //!< cutoff distance
  ActiveNematicParameters** m_pair_params;   //!< type specific pair parameters 
    
};

typedef shared_ptr<PairActiveNematicPotential> PairActiveNematicPotentialPtr;

#endif
