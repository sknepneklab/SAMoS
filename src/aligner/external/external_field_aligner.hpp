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
 *   (c) 2013, 2014, 2015
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file external_field_aligner.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 11-Apr-2015
 * \brief Declaration of ExternalFieldAlign class
 */ 

#ifndef __EXTERNAL_FIELD_ALIGN_HPP__
#define __EXTERNAL_FIELD_ALIGN_HPP__

#include <cmath>
#include <vector>

#include "external_aligner.hpp"

using std::make_pair;
using std::sqrt;
using std::vector;

//! Structure that handles parameters for the polar alignment 
struct FieldAlignParameters
{
  double hx, hy, hz;
};


/*! ExternalFieldAlign implements the alignment to the external field 
 *  given as constant vector \f$ \vec H \f$. Torque on particle \f$ i \f$ is then 
 *  \f$ \vec \tau_i = \vec H \times \vec n_i \f$. 
 */
class ExternalFieldAlign : public ExternalAlign
{
public:
  
   //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters 
  ExternalFieldAlign(SystemPtr sys, MessengerPtr msg, pairs_type& param) : ExternalAlign(sys,msg,param)
  {
    int ntypes = m_system->get_ntypes();
    if (param.find("hx") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External field aligner. No x component of external field given. Assuming 0.");
      m_hx = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External field aligner. x component of the external field set to "+param["hx"]+".");
      m_hx = lexical_cast<double>(param["hx"]);
    }
    if (param.find("hy") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External field aligner. No y component of external field given. Assuming 0.");
      m_hy = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External field aligner. y component of the external field set to "+param["hy"]+".");
      m_hy = lexical_cast<double>(param["hy"]);
    }
    if (param.find("hz") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External field aligner. No z component of external field given. Assuming 0.");
      m_hz = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External field aligner. z component of the external field set to "+param["hz"]+".");
      m_hz = lexical_cast<double>(param["hz"]);
    }
    
    m_type_params = new FieldAlignParameters[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
        m_type_params[i].hx = m_hx;
        m_type_params[i].hy = m_hy;
        m_type_params[i].hz = m_hz;
    }
  }
  
  virtual ~ExternalFieldAlign()
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
      m_msg->msg(Messenger::ERROR,"Type has not been defined for external field alignment.");
      throw runtime_error("Missing key for external alignment parameters.");
    }
    
    type = lexical_cast<int>(self_param["type"]);
    
    
    if (self_param.find("hx") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External field aligner. Setting x component of external filed to "+self_param["hx"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["hx"] = lexical_cast<double>(self_param["hx"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External field aligner. Using default x component of external filed ("+lexical_cast<string>(m_hx)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["hx"] = m_hx;
    }
    if (self_param.find("hy") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External field aligner. Setting y component of external filed to "+self_param["hy"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["hy"] = lexical_cast<double>(self_param["hy"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External field aligner. Using default y component of external filed ("+lexical_cast<string>(m_hy)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["hy"] = m_hy;
    }
    if (self_param.find("hz") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External field aligner. Setting z component of external filed to "+self_param["hz"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["hz"] = lexical_cast<double>(self_param["hz"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External field aligner. Using default z component of external filed ("+lexical_cast<string>(m_hz)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["hz"] = m_hz;
    }
    
    m_type_params[type].hx = param["hx"];
    m_type_params[type].hy = param["hy"];
    m_type_params[type].hz = param["hz"];
    
    m_has_params = true;
  }
 
  
  //! Computes "torques"
  void compute();
  
  
private:
       
  double m_hx;                             //!< x component of the external field
  double m_hy;                             //!< y component of the external field
  double m_hz;                             //!< z component of the external field
  FieldAlignParameters* m_type_params;   //!< type specific pair parameters 
     
};

typedef shared_ptr<ExternalFieldAlign> ExternalFieldAlignPtr;

#endif

