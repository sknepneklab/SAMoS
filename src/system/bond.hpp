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
 * \file bond.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Nov-2014
 * \brief Declaration of Bond data structure.
 */ 

#ifndef __BOND_HPP__
#define __BOND_HPP__

const int NUM_BOND_ATTRIB = 4;  //!< Number of bond attributes (columns in the bond file)

/*! Bond structure
 *  This is a basic data structure that hold information about bonds
 *  in the system. This is used to model chain-like molecules, such as
 *  polymers.
 *  Bond holds two indices, i and j which are id's of particles
 *  that are bonded. It also includes bond type, which is a positive
 *  integer.
 */
struct Bond
{
  //! Construct a bond object
  //! \param id Unique id of the bond
  //! \param type Bond type
  //! \param i Id of the first particle in the pair
  //! \param j Id of the second particle in the pair
  Bond(int id, int type, int i, int j) : id(id), type(type), i(i), j(j)  { }
  int id;     //!< Unique id of a bond
  int type;   //!< Bond type
  int i;      //!< Id of the first particle in the pair
  int j;      //!< Id of the second particle in the pair
};


#endif