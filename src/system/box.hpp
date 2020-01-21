/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

/*!
 * \file box.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Simulation box
 */ 

#ifndef __BOX_H__
#define __BOX_H__

#include <cassert>
#include <memory>

using std::shared_ptr;

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
