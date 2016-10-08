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
 * \file edge.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 26-Nov-2015
 * \brief Auxiliary functions for Edge class
 */ 

#include "edge.hpp"

ostream& operator<<(ostream& out, const Edge& e)
{
  out << " ---------- EDGE ------------------ " << endl;
  out << format("id : %d\n") % e.id
  << format("from : %d\n") % e.from
  << format("to : %d\n") % e.to
  << format("face : %d\n") % e.face
  << format("pair : %d\n") % e.pair;
  if (e.boundary)
    out << "boundary edge" << endl;
  out << " ------------------------------------";
  out << endl;
  return out;
}

