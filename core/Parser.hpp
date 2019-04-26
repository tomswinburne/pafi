#ifndef PARSER_H
#define PARSER_H

#include <vector>
#include <map>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
#include <boost/foreach.hpp>



class Parser {
public:

// Loads config file, lammps scripts etc
Parser(std::string file);

// remove trailing whitespace - avoids RapidXML parsing bug
std::string rtws(std::string s);

void welcome_message();

std::vector<std::string> Parse(std::string r,bool needseed=true);

std::vector<std::string> Script(std::string sn);

void seed(int random_seed);

// Create empty property tree object
boost::property_tree::ptree tree;
std::map<std::string,std::string> parameters;
std::map<std::string,std::string> scripts;

std::vector<std::string> KnotList;
double lowT,highT,Friction;
int CoresPerWorker, nPlanes, TSteps;
std::string dump_dir;
bool seeded;

private:
  boost::random::mt11213b rng;
};

#endif // PARSER_H
