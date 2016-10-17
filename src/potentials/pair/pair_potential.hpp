/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file pair_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Declaration of PairPotential class
 */ 

#ifndef __PAIR_POTENTIAL_HPP__
#define __PAIR_POTENTIAL_HPP__

#include <map>
#include <vector>
#include <string>

#include "system.hpp"
#include "neighbour_list.hpp"
#include "value.hpp"
#include "parse_parameters.hpp"

using std::map;
using std::pair;
using std::string;
using std::vector;


/*! PairPotential is the base class (abstract) that handles
 *  all calls to the pair potential evaluations. Children 
 *  of this class will implement actual force fields.
 */
class PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in). Line 340 of samos.cpp: The m_val(val) is set with the number of time steps to phase in
  //! \param param Contains information about all parameters 
  PairPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : m_system(sys), 
                                                                                                            m_msg(msg),
                                                                                                            m_nlist(nlist),
                                                                                                            m_val(val),
                                                                                                            m_has_pair_params(false),
                                                                                                            m_shifted(false),
                                                                                                            m_use_particle_radii(false),
                                                                                                            m_phase_in(false),
                                                                                                            m_compute_stress(false)
  {  
    m_known_params.push_back("min_val");
    m_known_params.push_back("max_val");
    m_known_params.push_back("ntypes");
    if (param.find("min_val") != param.end())
    {
      if (lexical_cast<double>(param["min_val"]) != 0.0)
        m_msg->msg(Messenger::WARNING,"Minimum value for the phase-in factor is not set to 0 (it is equal "+param["min_val"]+"). Is that what you really want to do?");
    }
    if (param.find("max_val") != param.end())
    {
      if (lexical_cast<double>(param["max_val"]) != 1.0)
        m_msg->msg(Messenger::WARNING,"Maximum value for the phase-in factor is not set to 1 (it is equal "+param["max_val"]+"). Is that what you really want to do?");
    }
    if (param.find("ntypes") != param.end())
    {
      m_ntypes = lexical_cast<int>(param["ntypes"]);
      if (m_ntypes < m_system->get_ntypes())
      {
        m_msg->msg(Messenger::WARNING,"Number of particle types has to be equal or greater than the number of types defined in the intial configuation. Using value read from the intial configuration.");
        m_ntypes = m_system->get_ntypes();
      }
      else
        m_msg->msg(Messenger::INFO,"Number of different particle types is set to "+param["ntypes"]+".");
    }
    else
    {
      m_ntypes = m_system->get_ntypes();
      m_msg->msg(Messenger::WARNING,"Number of particle types will be read from the intitial configuration.");
    }
  }
                                                                                                       
  //! Destructor 
  virtual ~PairPotential() { }
                                                                                                   
  //! Get the total potential energy
  double get_potential_energy() { return m_potential_energy; } //!< \return value of the total potential energy
  
  //! Set pair parameters data for pairwise interactions    
  virtual void set_pair_parameters(pairs_type&) = 0;
  
  //! Set type parameters for each particle
  virtual void set_type_parameters(pairs_type&) { }  
  
  //! Returns true if the specific pair potential needs neighbour list
  virtual bool need_nlist() = 0;
  
  //! Computes potentials and forces for all particles
  virtual void compute(double) = 0;

  //! Check if there are no illegal parameters
  string params_ok(pairs_type& params)
  {
    for (pairs_type::iterator it_p = params.begin(); it_p != params.end(); it_p++)
      if (find(m_known_params.begin(),m_known_params.end(),(*it_p).first) == m_known_params.end())
        return (*it_p).first;
    return "";
  }
  
protected:
       
  SystemPtr m_system;               //!< Pointer to the System object
  MessengerPtr m_msg;               //!< Handles messages sent to output
  NeighbourListPtr m_nlist;         //!< Handles NeighbourList object
  ValuePtr m_val;                   //!< Value object for phasing in particles
  bool m_has_pair_params;           //!< Flag that controls if pair parameters are set
  bool m_shifted;                   //!< If true, potential is shifted at cutoff
  double m_potential_energy;        //!< Total potential energy
  bool m_use_particle_radii;        //!< If true, base potential ranges (if they exist) on particle radii
  bool m_phase_in;                  //!< If true, gradually switch on potential for particles that are younger than a given age
  bool m_compute_stress;            //!< If true, compute stress tensor
  int m_ntypes;                     //!< Total number of particle types in the system
  vector<string> m_known_params;    //!< Lists all known parameters accepted by a given pair potential
  
};


typedef shared_ptr<PairPotential> PairPotentialPtr;


#endif
