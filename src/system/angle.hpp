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
 * \file angle.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Nov-2014
 * \brief Declaration of Angle data structure.
 */ 

#ifndef __ANGLE_HPP__
#define __ANGLE_HPP__


const int NUM_ANGLE_ATTRIB = 5;  //!< Number of angle attributes (columns in the angle file)

/*! Angle structure
 *  This is a basic data structure that hold information about angle
 *  in the system. This is used to model bending rigidity 
 *  in chain-like molecules, such as polymers.
 *  Angle holds three indices, i, j and k which are id's of particles
 *  that are bonded. It also includes bond type, which is a positive
 *  integer.
 */
struct Angle
{
  //! Construct a angle object
  //! \param id Unique id of the bond
  //! \param i Id of the first (left) particle 
  //! \param j Id of the second (middle) particle 
  //! \param k If of the third (right) particle
  //! \param type Bond type
  Angle(int id, int type, int i, int j, int k) : id(id), type(type), i(i), j(j), k(k) { }
  int id;     //!< Unique id of a bond
  int type;   //!< Angle type
  int i;      //!< Id of the first (left) particle
  int j;      //!< Id of the second (middle) particle
  int k;      //!< If of the third (right) particle 
};


#endif