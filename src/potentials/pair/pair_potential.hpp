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
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters 
  PairPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : m_system(sys), 
                                                                                                            m_msg(msg),
                                                                                                            m_nlist(nlist),
                                                                                                            m_val(val),
                                                                                                            m_has_pair_params(false),
                                                                                                            m_shifted(false),
                                                                                                            m_use_particle_radii(false),
                                                                                                            m_phase_in(false)
                                                                                                            { }
                                                                                                       
  //! Destructor 
  virtual ~PairPotential() { }
  
                                                                                                   
  //! Get the total potential energy
  double get_potential_energy() { return m_potential_energy; } //!< \return value of the total potential energy
  
  //! Set pair parameters data for pairwise interactions    
  virtual void set_pair_parameters(pairs_type&) = 0;
  
  //! Returns true if the specific pair potential needs neighbour list
  virtual bool need_nlist() = 0;
  
  //! Computes potentials and forces for all particles
  virtual void compute(double) = 0;
  
protected:
       
  SystemPtr m_system;              //!< Pointer to the System object
  MessengerPtr m_msg;              //!< Handles messages sent to output
  NeighbourListPtr m_nlist;        //!< Handles NeighbourList object
  ValuePtr m_val;                  //!< Value object for phasing in particles
  bool m_has_pair_params;          //!< Flag that controls if pair parameters are set
  bool m_shifted;                  //!< If true, potential is shifted at cutoff
  double m_potential_energy;       //!< Total potential energy
  bool m_use_particle_radii;       //!< If true, base potential ranges (if they exist) on particle radii
  bool m_phase_in;                 //!< If true, gradually switch on potential for particles that are younger than a given age
  
};


typedef shared_ptr<PairPotential> PairPotentialPtr;


#endif
