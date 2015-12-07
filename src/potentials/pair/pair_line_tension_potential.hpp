/* *************************************************************
 *  
 *   VertexParticle Active Mater on Surfaces (SAMoS)
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
 * \file pair_line_tension_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 07-Dec-2015
 * \brief Declaration of PairLineTensionPotential class
 */ 

#ifndef __PAIR_LINE_TENSION_POTENTIAL_HPP__
#define __PAIR_LINE_TENSION_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;

//! Structure that handles parameters for the line tension pair potential
struct LineTensionParameters
{
  double lambda;
};

/*! PairLineTensionPotential implements the line tension term in the model for an active tissue model.
 *  Force only exists along boundary edges and is proportional to the lenght of the edge, i.e., 
 *  \f$ \vec f_i = \lambda \left|\vec r_{ij}\right|\hat\vec r_{ij} \f$ and \f$ \vec f_j -= \lambda \left|\vec r_{ij}\right|\hat\vec r_{ij} \f$,
 *  where \f$ \vec r_{ij} \f$ if the vector along the edge connecting vertices beloing to it.
 */
class PairLineTensionPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (k)
  PairLineTensionPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    int ntypes = m_system->get_ntypes();
    if (param.find("lambda") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No line tension (lambda) specified for line tension pair potential. Setting it to 1.");
      m_lambda = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global line tension (lambda) for line tension pair potential is set to "+param["lambda"]+".");
      m_lambda = lexical_cast<double>(param["lambda"]);
    }
    m_msg->write_config("potential.pair.line_tension.lambda",lexical_cast<string>(m_lambda));
    
    
    m_pair_params = new LineTensionParameters*[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
      m_pair_params[i] = new LineTensionParameters[ntypes];
      for (int j = 0; j < ntypes; j++)
      {
        m_pair_params[i][j].lambda = m_lambda;
      }
    }
    
  }
  
  virtual ~PairLineTensionPotential()
  {
    for (int i = 0; i < m_system->get_ntypes(); i++)
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
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in line tension potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in line tension potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    if (pair_param.find("lambda") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Line tension pair potential. Setting line tesion to "+pair_param["lambda"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["lambda"] = lexical_cast<double>(pair_param["lambda"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Line tension pair potential. Using default line tension ("+lexical_cast<string>(m_lambda)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["lambda"] = m_lambda;
    }
    m_msg->write_config("potential.pair.line_tension.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".push",lexical_cast<string>(param["lambda"]));
        
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
       
  double m_lambda;                  //!< line tension 
  LineTensionParameters** m_pair_params;       //!< type specific pair parameters 
     
};

typedef shared_ptr<PairLineTensionPotential> PairLineTensionPotentialPtr;

#endif
