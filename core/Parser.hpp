#ifndef PARSER_H
#define PARSER_H

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <random>

#include "rapidxml.hpp"


class Parser {
public:

// Loads config file, lammps scripts etc
Parser(std::string file);

// remove trailing whitespace - avoids RapidXML parsing bug
std::string rtws(std::string s);

std::string welcome_message();

std::vector<std::string> split_lines(std::string r);

std::vector<std::string> Script(std::string sn);

void seed(int random_seed);

std::string seed_str();

rapidxml::xml_document<> xml_doc;
rapidxml::xml_node<> * root_node;
rapidxml::xml_node<> * child_node;

std::map<std::string,std::string> parameters;
std::map<std::string,std::string> scripts;

std::vector<std::string> KnotList;
double lowT,highT,Friction,startr,stopr,maxjump_thresh,redo_thresh;
int CoresPerWorker, nPlanes, TSteps, nRepeats, maxExtraRepeats;
std::string dump_dir;
bool seeded,loglammps,postDump,preMin,xml_success;

private:
  std::mt19937 rng;
};

#endif // PARSER_H
