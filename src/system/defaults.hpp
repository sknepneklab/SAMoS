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
 * \file defaults.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 05-Nov-2013
 * \brief Some default constants (just to keep the all in one place)
 */ 

#ifndef __DEFAULTS_H__
#define __DEFAULTS_H__

#define DEFAULT_MESSENGER "messenger.msg"
#define DEFAULT_CUTOFF 3.0
#define DEFAULT_PADDING 0.5
#define DEFAULT_LX 10.0
#define DEFAULT_LY 10.0
#define DEFAULT_LZ 10.0
#define SMALL_NUMBER 1e-12
#define MAX_FACE 10
#ifndef NDEBUG
#define PRINT_EVERY 1
#else
#define PRINT_EVERY 10000
#endif

#endif
