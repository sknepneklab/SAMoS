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
 *   (c) 2013, 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file pair_nematic_aligner.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 05-Mar-2014
 * \brief Declaration of PairNemticAlign class
 */ 

#ifndef __PAIR_NEMATIC_ALIGN_HPP__
#define __PAIR_NEMATIC_ALIGN_HPP__

#include <cmath>
#include <vector>

#include "pair_aligner.hpp"

using std::make_pair;
using std::sqrt;
using std::vector;

/*! PairNematicAlign implements the "nematic" type alignment between neighbouring particles.
 *  For all particles within the cutoff distance \f$ r_{cut} \left(\leq r_{nl}\right) \f$
 *  (\f$ r_{nl} \f$ being neighbour list cutoff distance) we compute torque on the particle as
 *  \f$ \vec \tau_i = 2 J \sum_j \left(\vec n_i\cdot\vec n_j\right)\vec n_i\times\vec n_j \f$,
 *  where \f$ J \f$ is the coupling constant and \f$ \vec n_j \f$ is director of j-th neighbour. 
 */
class PairNematicAlign : public PairAlign
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param param Contains information about all parameters (J and cutoff distance)
  PairNematicAlign(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, pairs_type& param) : PairAlign(sys, msg, nlist, param)
  {
    if (param.find("J") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No coupling constant (J) specified for MF alignment. Setting it to 1.");
      m_J = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global coupling constant (J) for MF alignment is set to "+param["J"]+".");
      m_J = lexical_cast<double>(param["J"]);
    }
    if (param.find("rcut") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No cutoff distance (rcut) specified for MF alignment. Setting it to the global neighbour list cutoff.");
      m_rcut = m_nlist->get_cutoff();
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global cutoff distance (rcut) for MF alignment is set to "+param["rcut"]+".");
      m_rcut = lexical_cast<double>(param["rcut"]);
    }
  }
                                                                                                                
  //! Set pair parameters data for pairwise alignment    
  void set_pair_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    int type_1, type_2;
    
    if (pair_param.find("type_1") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pairwise alignment in MF aligner.");
      throw runtime_error("Missing key for pair alignment parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pairwise alignment in MF aligner.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    if (pair_param.find("J") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"MF pairwise alignment. Setting coupling constant to "+pair_param["J"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["J"] = lexical_cast<double>(pair_param["J"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"MF pairwise alignment. Using default strength ("+lexical_cast<string>(m_J)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["J"] = m_J;
    }
    if (pair_param.find("rcut") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"MF pairwise alignment. Setting cutoff distance to "+pair_param["rcut"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = lexical_cast<double>(pair_param["rcut"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"MF pairwise alignment. Using default cutoff distance ("+lexical_cast<string>(m_rcut)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = m_rcut;
    }    
        
    m_pair_params[make_pair(type_1,type_2)]["J"] = param["J"];
    if (type_1 != type_2)
      m_pair_params[make_pair(type_2,type_1)]["J"] = param["J"];
    
    m_pair_params[make_pair(type_1,type_2)]["rcut"] = param["rcut"];
    if (type_1 != type_2)
      m_pair_params[make_pair(type_2,type_1)]["rcut"] = param["rcut"];
    
    m_has_pair_params = true;
  }
  
  //! Returns true since MF alignment needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes "torques"
  void compute();
  
  
private:
       
  double m_J;       //!< Coupling constant
  double m_rcut;    //!< Cutoff distance (has to be less than neighbour list cutoff)
     
};

typedef shared_ptr<PairNematicAlign> PairNematicAlignPtr;

#endif
