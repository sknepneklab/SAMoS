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
 * \file box.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Simulation box
 */ 

#ifndef __BOX_H__
#define __BOX_H__

#include <cassert>
#include <boost/shared_ptr.hpp>

using boost::shared_ptr;

//! Simulation box data 
/*! Defines a rectangular simulation box. 
 *  For closed surfaces this is completely ignored.
 *  For systems in the plane this is important to enforce 
 *  proper boundary conditions.
 * 
*/
struct Box
{
  //! Construct the box
  /*!
	\param XLO x-coordinate of the low corner
	\param XHI x-coordinate of the high corner
	\param YLO y-coordinate of the low corner
	\param YHI y-coordinate of the high corner
	\param ZLO z-coordinate of the low corner
	\param ZHI z-coordinate of the high corner	
  */
  Box(double lx, double ly, double lz) : Lx(lx), Ly(ly), Lz(lz)
  {
    assert(Lx >= 0.0 && Ly >= 0.0 && Lz >= 0.0);
    xlo = -0.5*Lx;  xhi = 0.5*Lx;
    ylo = -0.5*Ly;  yhi = 0.5*Ly;
    zlo = -0.5*Lz;  zhi = 0.5*Lz;
  }
  Box(double XLO, double XHI, double YLO, double YHI, double ZLO, double ZHI) : xlo(XLO), xhi(XHI), ylo(YLO), yhi(YHI), zlo(ZLO), zhi(ZHI) 
  {
    assert(XHI > XLO && YHI > YLO && ZHI > ZLO);
    Lx = xhi - xlo;
    Ly = yhi - ylo;
    Lz = zhi - zlo;
  }
  //! Rescale Box
  //! \param a amount by which to rescale box
  void rescale(double a)
  {
    Lx *= a; Ly *= a; Lz *= a;
    xlo *= a;  ylo *= a;  zlo *= a;
    xhi *= a;  yhi *= a;  zhi *= a;
  }
  
  double xlo;         //!< low value of the x-coordinates of the box
  double xhi;         //!< high value of the x-coordinates of the box
  double ylo;         //!< low value of the y-coordinates of the box
  double yhi;         //!< high value of the y-coordinates of the box
  double zlo;         //!< low value of the z-coordinates of the box
  double zhi;         //!< high value of the z-coordinates of the box
  double Lx;          //!< box length in x direction
  double Ly;          //!< box length in y direction
  double Lz;          //!< box length in z direction
};

typedef shared_ptr<Box> BoxPtr;

#endif