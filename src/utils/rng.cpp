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
 * \file rng.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 24-Oct-2013
 * \brief Wrappers for the GSL random number generate
 */ 

#include "rng.hpp"

//! Initialize RNG (GSL)
//! \param seed initial seed for the random number generator
RNG::RNG(int seed) : m_seed(seed)
{
  gsl_rng_env_setup();
  GSL_RANDOM_TYPE = gsl_rng_default;
  GSL_RANDOM_GENERATOR = gsl_rng_alloc(GSL_RANDOM_TYPE);
  gsl_rng_set(GSL_RANDOM_GENERATOR, static_cast<unsigned long int>(m_seed));
}

//! Free data structures used for the random number generator
//! In case of GLS library we need to free some memory
RNG::~RNG()
{
  gsl_rng_free(GSL_RANDOM_GENERATOR);
}

//! Get a random number between 0 and 1 drawn from an uniform distribution
//! \return random number between 0 and 1
double RNG::drnd()
{
  return gsl_rng_uniform(GSL_RANDOM_GENERATOR);
}

//! Return a random number from a Gaussian distribution with a given standard deviation 
//! \param sigma standard deviation 
double RNG::gauss_rng(double sigma)
{
  return gsl_ran_gaussian(GSL_RANDOM_GENERATOR, sigma);
}

//! Get an integer random number between 0 and N drawn from an uniform distribution
//! \param N upper bound for the interval
//! \return integer random number between 0 and N
int RNG::lrnd(int N)
{
  return static_cast<int>(N*drnd());
}



