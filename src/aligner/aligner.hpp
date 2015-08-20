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
 * \file aligner.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Jan-2014
 * \brief Declaration of Aligner class.
 */ 

#ifndef __ALIGNER_HPP__
#define __ALIGNER_HPP__

#include <map>
#include <string>

using std::map;
using std::string;

#include "messenger.hpp"
#include "system.hpp"

#include "pair_aligner.hpp"
#include "external_aligner.hpp"

typedef map<string,PairAlignPtr> PairAlignType;
typedef map<string,ExternalAlignPtr> ExternAlignType;

/*! Aligner class handles all director aligners (pairwise and external field type) present 
 *  in the system. All aligners are stored in two STL maps which both have strings as keys
 *  while items are pointers (boost shared_ptr) to the PairAlign and ExternalAlign abstract classes,
 *  respectively. The integrator part of the code calls member functions of this class
 *  to perform actual alignment computation.
*/
class Aligner
{
public:
  
  //! Construct Aligner object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  Aligner(SystemPtr sys, const MessengerPtr msg) : m_system(sys), m_msg(msg), m_need_nlist(false) { }
  
  //! Destructor
  ~Aligner()
  {
    m_pair_align.clear();
    m_external_align.clear();
  }
  
  //! Add pairwise alignment to the list of all pair alignments
  //! \param name Unique name of the aligner 
  //! \param align Pointer to the pairwise alignment object
  void add_pair_align(const string& name, PairAlignPtr align)
  {
    m_pair_align[name] = align;
    m_need_nlist = m_need_nlist | align->need_nlist();
    m_msg->msg(Messenger::INFO,"Added pairwise alignment : " + name + " to the list of pair alignments.");
    if (align->need_nlist())
      m_msg->msg(Messenger::INFO,"Pairwise alignment " + name + " has neighbour list. Neighbour list updates will be performed during the simulation.");
  }
  
  //! Add external alignment to the list of all external alignments 
  //! \param name Unique name of the aligner
  //! \param align Pointer to the external alignment object
  void add_external_align(const string& name, ExternalAlignPtr align)
  {
    m_external_align[name] = align;
    m_msg->msg(Messenger::INFO,"Added external aligner : " + name + " to the list of external alignments.");
  }
  
  //! Add pairwise alignment parameters
  //! \param name Unique name of the aligner 
  //! \param params maps with new aligner parameters
  void add_pair_align_parameters(const string& name, pairs_type& params)
  {
    m_pair_align[name]->set_pair_parameters(params);
  }
  
  //! Add external aligner parameters
  //! \param name Unique name of the aligner 
  //! \param params maps with new aligner parameters
  void add_external_align_parameters(const string& name, pairs_type& params)
  {
    m_external_align[name]->set_parameters(params);
  }
  
  //! Compute total pair alignment energy of a given type (for measurement)
  //! \param type pair alignment type
  double compute_pair_alignment_energy_of_type(const string& type)
  {
    if (m_pair_align.find(type) == m_pair_align.end())
    {
      m_msg->msg(Messenger::ERROR,"Trying to compute pair alignment of type " + type + " that is not defined for this system.");
      throw runtime_error("Pair alignment of type " + type + " not defined.");
    }
    return m_pair_align[type]->get_potential_energy();
  }
  
  //! Returns true if any of the potentials need neighbour list
  bool need_nlist() { return m_need_nlist; } 
  
  //! Compute all alignments
  void compute();
  
private:
  
  SystemPtr m_system;            //!< Contains pointer to the System object
  MessengerPtr m_msg;            //!< Handles messages sent to output
  
  PairAlignType m_pair_align;        //!< Contains information about all pair alignment 
  ExternAlignType m_external_align;  //!< Contains information about all external alignment
  
  bool m_need_nlist;                  //!< If true, there are potentials that need neighbour list
   
};

typedef shared_ptr<Aligner> AlignerPtr;

#endif