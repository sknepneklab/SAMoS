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
 * \file pair_vicsek_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 10-Dec-2013
 * \brief Declaration of PairVicsekPotential class
 */ 

#ifndef __PAIR_VICSEK_POTENTIAL_HPP__
#define __PAIR_VICSEK_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;

/*! PairVicsekPotential is not really a potential in the strict sense of the word
 *  but is instead an alignment rule that computes the average direction of the 
 *  velocity averaged over a certain neighbourhood of the particle.
 *  \f$ \langle \vec v_i\rangle = \sum_{\mathrm{j n.n i}} \vec v_j \f$ 
 */
class PairVicsekPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param param Contains information about all parameters (k)
  PairVicsekPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, pairs_type& param) : PairPotential(sys, msg, nlist, param)
  {
    if (param.find("v0") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No velocity magnitude (v0) specified for Vicsek pair potential. Setting it to 1.");
      m_v0 = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Velocity magnitude (v0) for Vicsek pair potential is set to "+param["v0"]+".");
      m_v0 = lexical_cast<double>(param["v0"]);
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
  
  //! Computes potentials and forces for all particles
  void compute();
  
  
private:
       
  double m_v0;       //!< velocity magnitude
     
};

typedef shared_ptr<PairVicsekPotential> PairVicsekPotentialPtr;

#endif
