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
 * \file rng.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 24-Oct-2013
 * \brief Class RNG provides wrappers for the GSL random number generate
 */ 

#ifndef __RNG_H__
#define __RNG_H__

#include <memory>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using std::shared_ptr;

/*! Class handles random numbers in the system */
class RNG
{
public:
  
  //! Constructor (initialize random number generator)
  RNG(int);
  
  //! Destructor
  ~RNG();
  
  //! Return random number between 0 and 1
  double drnd();
  
  //! Return random integer between 0 and N
  int lrnd(int);
  
  //! Return a Gaussian distributed number with a given standard deviation
  double gauss_rng(double);

private:
  
  int m_seed;   //!< Random number generator seed
  const gsl_rng_type* GSL_RANDOM_TYPE;  //!< Pointer to the gsl_rng_type structure which handles the RNG type
  gsl_rng* GSL_RANDOM_GENERATOR;        //!< Pointer which holds the actual random number generator

};

typedef shared_ptr<RNG> RNGPtr;

#endif
