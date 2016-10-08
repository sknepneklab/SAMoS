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
 * \file pair_boundary_bending_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-04-2016
 * \brief Declaration of PairBoundaryBendingPotential class
 */ 

#ifndef __PAIR_BOUNDARY_BENDING_POTENTIAL_HPP__
#define __PAIR_BOUNDARY_BENDING_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;
using std::acos;

//! Structure that handles parameters for the boundary bending pair potential
struct BoundaryBendingParameters
{
  double kappa;
  double theta0;
};

/*! PairBoundaryBendingPotential implements the bending penalty to bend the boundary of a cell patch.
 *  It depends on the mesh being present in the simulation.
 */
class PairBoundaryBendingPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (k)
  PairBoundaryBendingPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param), m_kappa(1.0), m_theta0(M_PI), m_has_part_params(false)
  {
    m_known_params.push_back("kappa");
    m_known_params.push_back("theta");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for boundary bending potential.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in boundary potential.");
    }
    if (param.find("kappa") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No bending rigidity (kappa) specified for boundary bending pair potential. Setting it to 1.");
      m_kappa = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global boundary bending rigidity (kappa) is set to "+param["kappa"]+".");
      m_kappa = lexical_cast<double>(param["kappa"]);
    }
    m_msg->write_config("potential.pair.boundary_bending.kappa",lexical_cast<string>(m_kappa));
    
    if (param.find("theta") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No preferred angle (theta0) specified for boundary bending pair potential. Setting it to PI.");
      m_theta0 = M_PI;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global boundary prefered angle (theta0) is set to "+param["theta"]+".");
      m_theta0 = lexical_cast<double>(param["theta"]);
    }
    m_msg->write_config("potential.pair.boundary_bending.theta0",lexical_cast<string>(m_theta0));
    
    
    m_particle_params = new BoundaryBendingParameters[m_ntypes];
    for (int i = 0; i < m_ntypes; i++)
    {
      m_particle_params[i].kappa = m_kappa;
      m_particle_params[i].theta0 = m_theta0;
    }
    
  }
  
  virtual ~PairBoundaryBendingPotential()
  {
    delete [] m_particle_params;
  }
                                                                                                                
  //! Set type parameters data for individual particles
  void set_type_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    int type;
    
    if (pair_param.find("type") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type has not been defined for pair potential particle specific parameters in boundary bending potential.");
      throw runtime_error("Missing key for pair potential particle parameters.");
    }
    
    type = lexical_cast<int>(pair_param["type"]);
    
    if (pair_param.find("kappa") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Boundary bending pair potential. Setting bending rigidity to "+pair_param["kappa"]+" for particle pair of type "+lexical_cast<string>(type)+".");
      param["kappa"] = lexical_cast<double>(pair_param["kappa"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Boundary bending pair potential. Using default line tension ("+lexical_cast<string>(m_kappa)+") for particle pair of types "+lexical_cast<string>(type)+".");
      param["kappa"] = m_kappa;
    }
    m_msg->write_config("potential.pair.boundary_bending.type_"+pair_param["type"]+".push",lexical_cast<string>(param["kappa"]));
        
    if (pair_param.find("theta") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Boundary bending pair potential. Setting preferd angle to "+pair_param["theta"]+" for particle pair of types "+lexical_cast<string>(type)+".");
      param["theta"] = lexical_cast<double>(pair_param["theta"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Boundary bending pair potential. Using default preferred angle ("+lexical_cast<string>(m_theta0)+") for particle pair of type "+lexical_cast<string>(type)+".");
      param["theta"] = m_theta0;
    }
    m_msg->write_config("potential.pair.boundary_bending.type_"+pair_param["type"]+".push",lexical_cast<string>(param["theta"]));
        
    m_particle_params[type-1].kappa = param["kappa"];
    m_particle_params[type-1].theta0 = param["theta"];
    
    m_has_part_params = true;
  }
  
  //! Set pair parameters 
  void set_pair_parameters(pairs_type& pair_param) { }
  
  //! Returns true since vertex-particle potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
       
  double m_kappa;                  //!< bending rigidity 
  double m_theta0;                 //!< native angle 
  bool m_has_part_params;          //!< true if type specific particle parameters are given
  BoundaryBendingParameters* m_particle_params;       //!< type specific particle parameters 
     
};

typedef shared_ptr<PairBoundaryBendingPotential> PairBoundaryBendingPotentialPtr;

#endif
