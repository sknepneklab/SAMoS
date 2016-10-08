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
 * \file parse_aux.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 21-Oct-2013
 * \brief Auxiliary structures for parsing command lines
 */ 

#ifndef __PARSE_AUX_HPP__
#define __PARSE_AUX_HPP__

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/mpl/print.hpp>
#include <boost/config/warning_disable.hpp>
#include <boost/detail/lightweight_test.hpp>
#include <boost/spirit/include/qi_char.hpp>
#include <boost/spirit/include/qi_string.hpp>
#include <boost/spirit/include/qi_nonterminal.hpp>
#include <boost/spirit/include/qi_numeric.hpp>
#include <boost/spirit/include/qi_action.hpp>
#include <boost/spirit/include/qi_operator.hpp>

#include <boost/spirit/repository/include/qi_distinct.hpp>
#include <boost/spirit/include/qi_no_skip.hpp>

#include <string>

/*! \note Adopted from http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp */

namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;

namespace spirit = boost::spirit;
namespace ascii = boost::spirit::ascii;
namespace repo = boost::spirit::repository;


//! Define metafunctions allowing to compute the type of the distinct()
//! and ascii::char_() constructs
namespace traits
{
  //! Metafunction allowing to get the type of any repository::distinct(...) 
  //! construct
  template <typename Tail>
  struct distinct_spec
  : spirit::result_of::terminal<repo::tag::distinct(Tail)>
  {};
  
  //! Metafunction allowing to get the type of any ascii::char_(...) construct
  template <typename String>
  struct char_spec
  : spirit::result_of::terminal<spirit::tag::ascii::char_(String)>
  {};
};

//! Define a helper function allowing to create a distinct() construct from 
//! an arbitrary tail parser
template <typename Tail>
inline typename traits::distinct_spec<Tail>::type
distinct_spec(Tail const& tail)
{
  return repo::distinct(tail);
}

//! Define a helper function allowing to create a ascii::char_() construct 
//! from an arbitrary string representation
template <typename String>
inline typename traits::char_spec<String>::type
char_spec(String const& str)
{
  return ascii::char_(str);
}

//! The following constructs the type of a distinct_spec holding a
//! charset("0-9a-zA-Z_") as its tail parser
typedef traits::char_spec<std::string>::type charset_tag_type;
typedef traits::distinct_spec<charset_tag_type>::type keyword_tag_type;

//! Define a new Qi 'keyword' directive usable as a shortcut for a
//! repository::distinct(char_(std::string("0-9a-zA-Z_")))
std::string const keyword_spec("0-9a-zA-Z_");
keyword_tag_type const keyword = distinct_spec(char_spec(keyword_spec)); 


#endif
