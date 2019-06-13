/* ***************************************************************************
 *
 *  Copyright (C) 2013-20119 University of Dundee
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
 * \file class_factory.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 13-Jun-2019
 * \brief Define class factory 
*/

#ifndef __CLASS_FACTORY_HPP__
#define __CLASS_FACTORY_HPP__

#include <memory>

using std::make_shared;

template<typename ptrT>
class factory
{
  public:

    factory() { }

    template <typename... Args>
    ptrT operator()(Args ... args) 
    {
      return make_shared<T>(args ...);
    }

  private:

    typedef typename ptrT::element_type T;

};

#endif
