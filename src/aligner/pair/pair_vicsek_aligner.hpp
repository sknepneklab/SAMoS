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
