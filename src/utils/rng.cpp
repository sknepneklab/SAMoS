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
 *   (c) 2013
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file rng.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 24-Oct-2013
 * \brief Wrappers for the GSL random number generate
 */ 

#include "rng.hpp"

//! Initialize RNG (GSL)
//! \param seed initial seed for the random number generator
RNG::RNG(int seed) : m_seed(seed);
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



