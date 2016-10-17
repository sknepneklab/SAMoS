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
 * \file parse_parameters.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 21-Oct-2013
 * \brief Grammar for parsing parameters 
 */ 

#ifndef __PARSE_PARAMETERS__
#define __PARSE_PARAMETERS__

#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <iostream>
#include <map>

/*! \note Adopted from: http://svn.boost.org/svn/boost/trunk/libs/spirit/example/qi/key_value_sequence.cpp */

namespace qi = boost::spirit::qi;

typedef std::map<std::string, std::string> pairs_type;  //!< Holds (key,value) pair

struct key_value_sequence : qi::grammar<std::string::iterator, pairs_type(), qi::space_type>
{
  key_value_sequence() : key_value_sequence::base_type(query)
  {
    query =  (qi::lit('{') >> qi::lit('}')) | (qi::lit('{') >> pair >> *((qi::lit(';') | qi::lit(',')) >> pair) >> *(qi::lit(';')) >> qi::lit('}'));
    pair  =  key >> -('=' >> value);
    key   =  qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9");
    value =  +qi::char_("-a-zA-Z_0-9+./");
    //value =  qi::char_("-0-9.") >> *qi::char_("0-9Ee+.-");
  }
  
  qi::rule<std::string::iterator, pairs_type(), qi::space_type> query;
  qi::rule<std::string::iterator, std::pair<std::string, std::string>(), qi::space_type> pair;
  qi::rule<std::string::iterator, std::string(), qi::space_type> key, value;
};


#endif
