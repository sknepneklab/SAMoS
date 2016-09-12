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
 * \file external_shape_aligner.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 12-Sep-2016
 * \brief Declaration of ExternalShapeAlign class
 */ 

#ifndef __EXTERNAL_SHAPE_ALIGN_HPP__
#define __EXTERNAL_SHAPE_ALIGN_HPP__

#include <cmath>
#include <vector>

#include "external_aligner.hpp"

using std::make_pair;
using std::sqrt;
using std::vector;

//! Structure that handles parameters for the shape alignment 
struct ShapeParameters
{
  double J;
};


/*! ExternalShapeAlign implements the single particle alignment to the "shape" of the cell
 *  in tissue simulations. By shape, we mean the direction of the eigenvector of the 
 *  cell's gyration tensor corresponding to the largest eigenvalue. 
 *  For all particles we compute torque on the particle as
 *  \f$ \vec \tau_i = -J \vec n_i\times\vec a_i \f$,
 *  where \f$ J \f$ is the coupling strength and \f$ \vec a_i \f$ is unit-lenght eigenvector of
 *  the gyration tensor corresponding to the largest eigenvalue. 
 */
class ExternalShapeAlign : public ExternalAlign
{
public:
  
   //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters 
  ExternalShapeAlign(SystemPtr sys, MessengerPtr msg, pairs_type& param) : ExternalAlign(sys,msg,param)
  {
    int ntypes = m_system->get_ntypes();
    if (param.find("J") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No coupling constant (J) specified for cell shape alignment. Setting it to 1.");
      m_J = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global coupling constant (J) for cell shape alignment is set to "+param["J"]+".");
      m_J = lexical_cast<double>(param["J"]);
    }
    m_msg->write_config("aligner.external.shape.J",lexical_cast<string>(m_J));
    
    m_type_params = new ShapeParameters[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
        m_type_params[i].J = m_J;
    }
  }
  
  virtual ~ExternalShapeAlign()
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
      m_msg->msg(Messenger::ERROR,"Type has not been defined for external alignment in cell shape aligner.");
      throw runtime_error("Missing key for external alignment parameters.");
    }
    
    type = lexical_cast<int>(self_param["type"]);
    
    if (self_param.find("J") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External cell shape aligner. Setting coupling constant to "+self_param["J"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["J"] = lexical_cast<double>(self_param["J"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External cell shape aligner. Using default strength ("+lexical_cast<string>(m_J)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["J"] = m_J;
    }
    m_msg->write_config("aligner.external.shape.type."+self_param["type"]+".J",lexical_cast<string>(param["J"]));
    
    m_type_params[type].J = param["J"];
        
    m_has_params = true;
  }
 
  
  //! Computes "torques"
  void compute();
  
  
private:
       
  double m_J;                             //!< Coupling constant
  ShapeParameters* m_type_params;        //!< type specific pair parameters 
     
};

typedef shared_ptr<ExternalShapeAlign> ExternalShapeAlignPtr;

#endif

