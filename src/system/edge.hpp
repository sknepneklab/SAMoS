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
 * \file edge.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Nov-2015
 * \brief Declaration of Edge class.
 */ 

#ifndef __EDGE_HPP__
#define __EDGE_HPP__

const int NO_FACE = -1;   // value indicting that one of the faces is not set

using std::endl;

/*! Edge class keeps track of the edge information in the mesh
 *
 */
struct Edge
{
  //! Construct a Edge object
  //! \param id edge id
  //! \param i id of 1st vertex
  //! \param j id of 2nd vertex
  Edge(int id, int i, int j) : id(id), i(i), j(j), f1(NO_FACE), f2(NO_FACE), boundary(false)  {   }
  
  
  int id;                      //!< Edge id
  int i;                       //!< 1st vertex
  int j;                       //!< 2nd vertex

  int f1;                      //!< 1st face (NO_FACE if none)
  int f2;                      //!< 2nd face (NO_FACE if none)
  
  bool boundary;               //!< If true, this is a boundary edge
    
};

#endif