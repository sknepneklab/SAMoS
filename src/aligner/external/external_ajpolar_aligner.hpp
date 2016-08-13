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
 * \file external_ajpolar_aligner.hpp
 * \author Silke Henkes, silkehenkes@gmail.com
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Oct-2014
 * \brief Declaration of ExternalAJPolarAlign class
 */ 

#ifndef __EXTERNAL_AJPOLAR_ALIGN_HPP__
#define __EXTERNAL_AJPOLAR_ALIGN_HPP__

#include <cmath>
#include <vector>

#include "external_aligner.hpp"

using std::make_pair;
using std::sqrt;
using std::vector;

//! Structure that handles parameters for the polar alignment 
struct PolarAJAlignParameters
{
  double tau;
};


/*! ExternalAJPolarAlign implements the active jamming type "polar" type alignment on a single particle.
 *  For all particles we compute torque on the particle as
 *  \f$ \vec \tau_i = -1/tau \vec n_i\times\vec v_i \f$,
 *  where \f$ tau \f$ is the coupling alignment time and \f$ \vec v_i \f$ is the velocity of the particle itself. 
 */
class ExternalAJPolarAlign : public ExternalAlign
{
public:
  
   //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters 
  ExternalAJPolarAlign(SystemPtr sys, MessengerPtr msg, pairs_type& param) : ExternalAlign(sys,msg,param)
  {
    int ntypes = m_system->get_ntypes();
    if (param.find("tau") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No coupling constant (tau) specified for active jamming polar alignment. Setting it to 1.");
      m_tau = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global coupling constant (tau) for active jamming polar alignment is set to "+param["tau"]+".");
      m_tau = lexical_cast<double>(param["tau"]);
    }
    
    if (param.find("normalise") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Velocity alignment is set to use normalised velocity vectors. Alignment strength is independent of v0.");
      m_msg->write_config("aligner.external.ajpolar.normalise","true");
      m_normalise = true;
    }
    else
    {
      m_msg->write_config("aligner.external.ajpolar.normalise","true");
      m_normalise = false;
    }
    m_msg->write_config("aligner.external.ajpolar.tau",lexical_cast<string>(m_tau));
    m_type_params = new PolarAJAlignParameters[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
        m_type_params[i].tau = m_tau;
    }
  }
  
  virtual ~ExternalAJPolarAlign()
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
      m_msg->msg(Messenger::ERROR,"Type has not been defined for external alignment in active jamming polar aligner.");
      throw runtime_error("Missing key for external alignment parameters.");
    }
    
    type = lexical_cast<int>(self_param["type"]);
    
    if (self_param.find("tau") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External active jamming polar aligner. Setting coupling constant to "+self_param["tau"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["tau"] = lexical_cast<double>(self_param["tau"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External active jamming polar aligner. Using default strength ("+lexical_cast<string>(m_tau)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["tau"] = m_tau;
    }
    m_msg->write_config("aligner.external.ajpolar.type."+self_param["type"]+".tau",lexical_cast<string>(param["tau"]));
    
    m_type_params[type].tau = param["tau"];
        
    m_has_params = true;
  }
 
  
  //! Computes "torques"
  void compute();
  
  
private:
       
  double m_tau;                             //!< Coupling constant
  bool m_normalise;                        //!< Whether to normalise the velocity in the alignment term
  PolarAJAlignParameters* m_type_params;   //!< type specific pair parameters 
     
};

typedef shared_ptr<ExternalAJPolarAlign> ExternalAJPolarAlignPtr;

#endif

