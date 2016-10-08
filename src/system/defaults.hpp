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
