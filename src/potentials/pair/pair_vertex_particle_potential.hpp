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
 * \file pair_vertex_particle_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-Nov-2015
 * \brief Declaration of PairVertexParticlePotential class
 */ 

#ifndef __PAIR_VERTEX_PARTICLE_POTENTIAL_HPP__
#define __PAIR_VERTEX_PARTICLE_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;

//! Structure that handles parameters for the vertex-particle pair potential
struct VertexParticleParameters
{
  double K;
  double gamma;
  double lambda;
};

/*! PairVertexParticlePotential implements the vertex-particle model for an active tissue model.
 *  \todo Document the force.
 */
class PairVertexParticlePotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (k)
  PairVertexParticlePotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param), m_has_part_params(false), m_include_boundary(false)
  {
    m_known_params.push_back("K");
    m_known_params.push_back("gamma");
    m_known_params.push_back("lambda");
    m_known_params.push_back("phase_in");
    m_known_params.push_back("compute_stress");
    m_known_params.push_back("mesh_update_steps");
    m_known_params.push_back("include_boundary");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for vertex-particle pair potential.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in vertex-particle model.");
    }
    if (param.find("K") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No area stiffness (K) specified for vertex-particle pair potential. Setting it to 1.");
      m_K = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global area stiffness (K) for vertex-particle pair potential is set to "+param["K"]+".");
      m_K = lexical_cast<double>(param["K"]);
    }
    m_msg->write_config("potential.pair.vertex_particle.K",lexical_cast<string>(m_K));
    if (param.find("gamma") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No perimiter stiffness (gamma) specified for vertex-particle pair potential. Setting it to 1.");
      m_gamma = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global perimeter stiffness (gamma) for vertex-particle pair potential is set to "+param["gamma"]+".");
      m_gamma = lexical_cast<double>(param["gamma"]);
    }
    m_msg->write_config("potential.pair.vertex_particle.gamma",lexical_cast<string>(m_gamma));
    if (param.find("lambda") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No contact stiffness (lambda) specified for vertex-particle pair potential. Setting it to 1.");
      m_lambda = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global contact stiffness (lambda) for vertex-particle pair potential is set to "+param["lambda"]+".");
      m_lambda = lexical_cast<double>(param["lambda"]);
    }
    m_msg->write_config("potential.pair.vertex_particle.lambda",lexical_cast<string>(m_lambda));
    
    if (param.find("phase_in") != param.end())
    {
      m_msg->msg(Messenger::INFO,"VertexParticle pair potential. Gradually phasing in the potential for new particles.");
      m_phase_in = true;
      m_msg->write_config("potential.pair.vertex_particle.phase_in","true");
    }    
    
    if (param.find("compute_stress") != param.end())
    {
      m_msg->msg(Messenger::INFO,"VertexParticle pair potential. Compute per cell stress tensor.");
      m_compute_stress = true;
      m_msg->write_config("potential.pair.vertex_particle.compute_stress","true");
    }
    
    if (param.find("mesh_update_steps") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Number of time steps between mesh updates in vertex-particle pair potential not set. Assuming default 0. No updates beyond those during the neighbour list and/or regular equiangulation steps for triangular lattices will be made.");
      m_mesh_update_steps = 0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Number of time steps between mesh updates in vertex-particle pair potential set to "+param["mesh_update_steps"]+".");
      m_mesh_update_steps = lexical_cast<int>(param["mesh_update_steps"]);
    }
    m_msg->write_config("potential.pair.vertex_particle.mesh_update_steps",lexical_cast<string>(m_mesh_update_steps));
    
    if (param.find("include_boundary") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"Vertex-particle pair potential INCLUDES additional terms arising at the boundary. This is a future feature, likely not to work.");
      m_include_boundary = true;
      m_msg->write_config("potential.pair.vertex_particle.include_boundary","true");
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Vertex-particle pair potential does not include additional terms arising at the boundary.");
      m_msg->write_config("potential.pair.vertex_particle.include_boundary","false");
    }
    
    m_particle_params = new VertexParticleParameters[m_ntypes];
    for (int i = 0; i < m_ntypes; i++)
    {
      m_particle_params[i].K = m_K;
      m_particle_params[i].gamma = m_gamma;
      m_particle_params[i].lambda = m_lambda;
    }
    
    m_pair_params = new VertexParticleParameters*[m_ntypes];
    for (int i = 0; i < m_ntypes; i++)
    {
      m_pair_params[i] = new VertexParticleParameters[m_ntypes];
      for (int j = 0; j < m_ntypes; j++)
      {
        m_pair_params[i][j].K = m_K;
        m_pair_params[i][j].gamma = m_gamma;
        m_pair_params[i][j].lambda = m_lambda;
      }
    }
    
  }
  
  virtual ~PairVertexParticlePotential()
  {
    for (int i = 0; i < m_ntypes; i++)
      delete [] m_pair_params[i];
    delete [] m_pair_params;
    delete [] m_particle_params;
  }
  
  //! Set type parameters data for individual particles
  void set_type_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    int type;
    
    if (pair_param.find("type") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type has not been defined for type specific parameters in vertex-particle potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    
    type = lexical_cast<int>(pair_param["type"]);
        
    if (pair_param.find("K") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"VertexParticle pair potential. Setting area stiffness to "+pair_param["K"]+" for particles of type "+lexical_cast<string>(type)+".");
      param["K"] = lexical_cast<double>(pair_param["K"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"VertexParticle pair potential. Using default area stiffness ("+lexical_cast<string>(m_K)+") for particles of type "+lexical_cast<string>(type)+".");
      param["K"] = m_K;
    }
    m_msg->write_config("potential.pair.vertex_particle.type_"+pair_param["type"]+".push",lexical_cast<string>(param["K"]));
    if (pair_param.find("gamma") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"VertexParticle pair potential. Setting perimeter stiffness to "+pair_param["gamma"]+" for particles of type "+lexical_cast<string>(type)+".");
      param["gamma"] = lexical_cast<double>(pair_param["gamma"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"VertexParticle pair potential. Using default perimeter stiffness ("+lexical_cast<string>(m_gamma)+") for particles of type "+lexical_cast<string>(type)+".");
      param["gamma"] = m_gamma;
    }
    m_msg->write_config("potential.pair.vertex_particle.type_"+pair_param["type"]+".push",lexical_cast<string>(param["gamma"]));
    if (pair_param.find("lambda") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"VertexParticle pair potential. Setting edge contractility (lambda) to "+pair_param["lambda"]+" for particles of type "+lexical_cast<string>(type)+".");
      param["lambda"] = lexical_cast<double>(pair_param["lambda"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"VertexParticle pair potential. Using default edge contractility ("+lexical_cast<string>(m_lambda)+") for particles of type "+lexical_cast<string>(type)+".");
      param["lambda"] = m_lambda;
    }
    m_msg->write_config("potential.pair.vertex_particle.type_"+pair_param["type"]+".push",lexical_cast<string>(param["lambda"]));
    
        
    m_particle_params[type-1].K = param["K"];
    m_particle_params[type-1].gamma = param["gamma"];
    m_particle_params[type-1].lambda = param["lambda"];
        
    m_has_part_params = true;
  }
                                                                                                                
  //! Set pair parameters data for pairwise interactions    
  void set_pair_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    int type_1, type_2;
    
    if (pair_param.find("type_1") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in vertex-particle potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in vertex-particle potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    if (pair_param.find("lambda") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"VertexParticle pair potential. Setting contact stiffness to "+pair_param["lambda"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["lambda"] = lexical_cast<double>(pair_param["lambda"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"VertexParticle pair potential. Using default contact stiffness ("+lexical_cast<string>(m_lambda)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["lambda"] = m_lambda;
    }
    m_msg->write_config("potential.pair.vertex_particle.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".push",lexical_cast<string>(param["lambda"]));
        
    m_pair_params[type_1-1][type_2-1].lambda = param["lambda"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].lambda = param["lambda"];
    
    m_has_pair_params = true;
  }
  
  //! Returns true since vertex-particle potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
       
  double m_K;                       //!< cell area stiffness
  double m_gamma;                   //!< cell perimeter stiffness
  double m_lambda;                  //!< cell contact stiffness
  bool m_has_part_params;           //!< true if type specific particle parameters are given
  int m_mesh_update_steps;          //!< number of time steps between mesh (tessalation) recalculation
  bool m_include_boundary;          //!< if true, include boudary terms in force calculation
  VertexParticleParameters*  m_particle_params;   //!< type specific particle parameters 
  VertexParticleParameters** m_pair_params;       //!< type specific pair parameters 
     
};

typedef shared_ptr<PairVertexParticlePotential> PairVertexParticlePotentialPtr;

#endif
