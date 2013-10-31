#include "parse_command.hpp"
#include "parse_parameters.hpp"
#include "parse_potential.hpp"
#include "parse_rng_seed.hpp"
#include "parse_input.hpp"

#include <iostream>
#include <string>
#include <map>

typedef std::map<std::string, qi::grammar<std::string::iterator, qi::space_type>* > parse_type;

int main(int argc, char* argv[])
{
  //std::string s = "pair_potential lj { eps = 1.0, sigma = 1.0 }";
  std::string s = "input meho-je.silni";
  parse_type command;
  CommandData command_data;
  PotentialData pot_data;
  RngSeedData rng_seed_data;
  InputData inp_data;
  pairs_type m;
 
  
  command_grammar cg(command_data);
  
  potential_grammar pg(pot_data);
  
  rng_seed_grammar rng(rng_seed_data);
  
  input_grammar inp(inp_data);
  
  key_value_sequence ks;
  
  command["pair_potential"] = &pg;
  command["seed"] = &rng;
  command["input"] = &inp;
  
  if (!qi::phrase_parse(s.begin(), s.end(), cg, qi::space))
  {
    std::cout << "-------------------------------- \n";
    std::cout << "Command parsing failed\n";
    std::cout << "-------------------------------- \n";
  }
  else
  {
    std::cout << command_data.command << std::endl;
    std::cout << command_data.attrib_param_complex << std::endl;
    
    if (!qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), *(command[command_data.command]) , qi::space))
    {
      std::cout << "-------------------------------- \n";
      std::cout << "Parsing failed for command " << command_data.command << std::endl;
      std::cout << "-------------------------------- \n";
    }
    else
    {
      if (command_data.command == "pair_potential")
      {
        std::cout << pot_data.params << std::endl;
        if (!qi::phrase_parse(pot_data.params.begin(), pot_data.params.end(), ks , qi::space, m))
        {
          std::cout << "-------------------------------- \n";
          std::cout << "Pair Potential parameters failed\n";
          std::cout << "-------------------------------- \n";
        }
        else
        {
          std::cout << "-------------------------------- \n";
          std::cout << "Parsing succeeded, found entries:\n";
          pairs_type::iterator end = m.end();
          for (pairs_type::iterator it = m.begin(); it != end; ++it)
          {
            std::cout << (*it).first;
            if (!(*it).second.empty())
              std::cout << "=" << (*it).second;
            std::cout << std::endl;
          }
          std::cout << "---------------------------------\n";
          
        }
      }
      else if (command_data.command == "seed")
      {
        std::cout << "-------------------------------- \n";
        std::cout << "Parsing succeeded, found entries:\n";
        std::cout << rng_seed_data.seed << std::endl;
      }
      else if (command_data.command == "input")
      {
        std::cout << "-------------------------------- \n";
        std::cout << "Parsing succeeded, found entries:\n";
        std::cout << inp_data.name << std::endl;
      }
    }
    
  }
  
  
  return 0;
}