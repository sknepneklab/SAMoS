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
    query =  qi::lit('{') >> pair >> *((qi::lit(';') | qi::lit(',')) >> pair) >> qi::lit('}');
    pair  =  key >> -('=' >> value);
    key   =  qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9");
    value =  qi::char_("-0-9.") >> *qi::char_("0-9Ee+.-");
  }
  
  qi::rule<std::string::iterator, pairs_type(), qi::space_type> query;
  qi::rule<std::string::iterator, std::pair<std::string, std::string>(), qi::space_type> pair;
  qi::rule<std::string::iterator, std::string(), qi::space_type> key, value;
};


#endif