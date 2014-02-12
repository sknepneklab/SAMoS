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
  m_system->reset_torques();
  PairAlignType::iterator it_pair;
  ExternAlignType::iterator it_ext;
  
  for(it_pair = m_pair_align.begin(); it_pair != m_pair_align.end(); it_pair++)
    (*it_pair).second->compute();
  for(it_ext = m_external_align.begin(); it_ext != m_external_align.end(); it_ext++)
    (*it_ext).second->compute();
}