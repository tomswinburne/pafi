#ifndef PARSER_H
#define PARSER_H

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>

#include "rapidxml.hpp"

std::string lammps_int_to_date(int stamp);

class Parser {
public:

// Loads config file, lammps scripts etc
Parser(std::string file, bool test);

// remove trailing whitespace - avoids RapidXML parsing bug
std::string rtws(std::string s);

std::string welcome_message();

void overwrite_xml(int nProcs);

void set_parameters();

std::vector<std::string> split_lines(std::string r);

std::vector<std::string> Script(std::string sn);

std::vector<double> sample_r(std::vector<double> pathway_r);

void seed(unsigned _random_seed);

std::string seed_str(bool reseed=true);

rapidxml::xml_document<> xml_doc;
rapidxml::xml_node<> * root_node;
rapidxml::xml_node<> * child_node;

std::map<std::string,std::string> parameters;
std::map<std::string,std::string> scripts;

std::vector<std::string> PathwayConfigurations;
std::vector<double> custom_positions;
double lowT,highT,Friction,startr,stopr,maxjump_thresh,redo_thresh,f_error_thresh;
int CoresPerWorker, nPlanes, TSteps, nRepeats, maxExtraRepeats, real_coord;
unsigned random_seed;
std::string dump_dir;
bool seeded,loglammps,postMin,preMin,xml_success;
bool workerDump,rediscretize,use_custom_positions,cubic_spline;

private:
  std::mt19937 rng;
};

#endif // PARSER_H
