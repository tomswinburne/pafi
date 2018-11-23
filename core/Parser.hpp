#ifndef PARSER_H
#define PARSER_H

#include <vector>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random.hpp>
#include <boost/foreach.hpp>


class Parser {
public:

// Loads config file, lammps scripts etc
Parser(std::string file);

// remove trailing whitespace - avoids RapidXML parsing bug
std::string rtws(std::string s);


std::vector<std::string> Parse(std::string r);

std::vector<std::string> Script(std::string sn);

// Create empty property tree object
boost::property_tree::ptree tree;
std::map<std::string,std::string> parameters;
std::map<std::string,std::string> scripts;

std::vector<std::string> KnotList;
int CoresPerWorker, nPlanes;

private:
  boost::random::mt11213b rng;
};

#endif // PARSER_H
