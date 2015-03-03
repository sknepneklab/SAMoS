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
 * \file regster_pair_aligners.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief Register all pair aligners with the class factory.
*/

#include "register.hpp"

void register_pair_aligners(PairAlignerMap& pair_aligners)
{
  // Register polar aligner with the pairwise aligner class factory
  pair_aligners["polar"] = boost::factory<PairPolarAlignPtr>();
  // Register nematic aligner with the pairwise aligner class factory
  pair_aligners["nematic"] = boost::factory<PairNematicAlignPtr>();
  // Register Vicsek aligner with the pairwise aligner class factory
  pair_aligners["vicsek"] = boost::factory<PairVicsekAlignPtr>();
}