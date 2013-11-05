/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
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

#include "parse_parameters.hpp"

using std::map;
using std::pair;
using std::string;
using std::vector;

typedef map<pair<int,int>, map<string,double> > PairData;

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
  //! \param param Contains information about all parameters 
  PairPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, pairs_type& param) : m_system(sys), 
                                                                                              m_msg(msg),
                                                                                              m_nlist(nlist),
                                                                                              m_has_pair_params(false),
                                                                                              m_shifted(false)
                                                                                              { }
                                                                                                                
  //! Get the total potential energy
  double get_potential_energy() { return m_potential_energy; } //!< \return value of the total potential energy
  
  //! Set pair parameters data for pairwise interactions    
  virtual void set_pair_parameters(pairs_type&) = 0;
  
  //! Computes potentials and forces for all particles
  virtual void compute() = 0;
  
  
protected:
       
  SystemPtr m_system;              //!< Pointer to the System object
  MessengerPtr m_msg;              //!< Handles messages sent to output
  NeighbourListPtr m_nlist;        //!< Handles NeighbourList object
  PairData m_pair_params;          //!< Handles specific parameters for a given pair
  bool m_has_pair_params;          //!< Flag that controls if pair parameters are set
  bool m_shifted;                  //!< If true, potential is shifted at cutoff
  double m_potential_energy;       //!< Total potential energy
  
};


typedef shared_ptr<PairPotential> PairPotentialPtr;


#endif
