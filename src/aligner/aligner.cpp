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
 * \file aligner.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Jan-2014
 * \brief Definition of Aligner class members.
*/ 

#include "aligner.hpp"

/*! Iterate over all pair and external aligners and compute 
 *  torques
 */
void Aligner::compute()
{
  //m_system->reset_torques();
  PairAlignType::iterator it_pair;
  ExternAlignType::iterator it_ext;
  
  for(it_pair = m_pair_align.begin(); it_pair != m_pair_align.end(); it_pair++)
    (*it_pair).second->compute();
  for(it_ext = m_external_align.begin(); it_ext != m_external_align.end(); it_ext++)
    (*it_ext).second->compute();
}