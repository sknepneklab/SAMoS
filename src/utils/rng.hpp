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
 * \file rng.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 24-Oct-2013
 * \brief Class RNG provides wrappers for the GSL random number generate
 */ 

#ifndef __RNG_H__
#define __RNG_H__

#include <boost/shared_ptr.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using boost::shared_ptr;

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
