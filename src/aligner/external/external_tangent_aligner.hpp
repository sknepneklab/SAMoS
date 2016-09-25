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
 *   (c) 2015, 2016
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015, 2016
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file external_tangent_aligner.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Sep-2016
 * \brief Declaration of ExternalTangentAlign class
 */ 

#ifndef __EXTERNAL_TANGENT_ALIGN_HPP__
#define __EXTERNAL_TANGENT_ALIGN_HPP__

#include <cmath>
#include <vector>

#include "external_aligner.hpp"

using std::make_pair;
using std::sqrt;
using std::vector;

//! Structure that handles parameters for the tangent to polymer alignment 
struct TangentAlignParameters
{
  double tau;
};


/*! ExternalTangentAlign implements the alignment of the director to the tangent to polymer.
 *  For a bead in the polymer chain we compute torque as 
 *  \f$ \vec \tau_i = -1/tau \sum_j\vec n_i\times\vec b_{ij} \f$,
 *  where \f$ tau \f$ is the coupling alignment time and \f$ \vec b_{ij} \f$ is the unit-length vector 
 *  pointing along the bond connecting partiles \f$ i \f$ and \f$ j $\f. 
 *  \note For this alignment to work, bonds have to be defined. User has to supply properly 
 *  ordered bonds in order to ensure that computed torques have desirable behaviour, e.g., to 
 *  align the director in the local direction of the polymer.   
 */
class ExternalTangentAlign : public ExternalAlign
{
public:
  
   //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters 
  ExternalTangentAlign(SystemPtr sys, MessengerPtr msg, pairs_type& param) : ExternalAlign(sys,msg,param)
  {
    if (m_system->num_bonds() == 0)
    {
      m_msg->msg(Messenger::ERROR,"Bonds have to defined in order to use aligner to the polymer tangent.");
      throw runtime_error("Bonds not defined in tangent aligner.");
    }
    int ntypes = m_system->get_ntypes();
    if (param.find("tau") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No coupling constant (tau) specified for tangent to polymer alignment. Setting it to 1.");
      m_tau = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global coupling constant (tau) for tangent to polymer alignment is set to "+param["tau"]+".");
      m_tau = lexical_cast<double>(param["tau"]);
    }
    m_msg->write_config("aligner.external.tangent.tau",lexical_cast<string>(m_tau));
    m_type_params = new TangentAlignParameters[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
        m_type_params[i].tau = m_tau;
    }
  }
  
  //! Destructor
  virtual ~ExternalTangentAlign()
  {
    delete [] m_type_params;
  }
      
   //! Set pair parameters data for self alignment   
  void set_parameters(pairs_type& self_param)
  {
    map<string,double> param;
    int type;
    
    if (self_param.find("type") == self_param.end())
    {
      m_msg->msg(Messenger::ERROR,"Type has not been defined for external alignment in tangent to polymer aligner.");
      throw runtime_error("Missing key for external alignment parameters.");
    }
    
    type = lexical_cast<int>(self_param["type"]);
    
    if (self_param.find("tau") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External tangent to polymer aligner. Setting coupling constant to "+self_param["tau"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["tau"] = lexical_cast<double>(self_param["tau"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External tangent to polymer aligner. Using default strength ("+lexical_cast<string>(m_tau)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["tau"] = m_tau;
    }
    m_msg->write_config("aligner.external.tangent.type."+self_param["type"]+".tau",lexical_cast<string>(param["tau"]));
    
    m_type_params[type].tau = param["tau"];
        
    m_has_params = true;
  }
 
  //! Computes "torques"
  void compute();
   
private:
       
  double m_tau;                             //!< Coupling constant
  bool m_normalise;                        //!< Whether to normalise the velocity in the alignment term
  TangentAlignParameters* m_type_params;   //!< type specific pair parameters 
     
};

typedef shared_ptr<ExternalTangentAlign> ExternalTangentAlignPtr;

#endif

