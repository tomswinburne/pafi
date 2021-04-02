#ifndef PARSER_H
#define PARSER_H

#include <stdio.h>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <random>
#include <tuple>
#include <iterator>
#include <regex>

#include "ConstantsTypes.hpp"
#include "rapidxml.hpp"

class Parser {
public:

// Loads config file, lammps scripts etc
Parser(std::string file, bool test);

// remove trailing whitespace - avoids RapidXML parsing bug
std::string rtws(std::string s);

std::string welcome_message();

void overwrite_xml(int nProcs);

void set_parameters();

void insert_params(std::string &s, Holder &params);

std::vector<std::string> split_lines(std::string r);

std::vector<std::string> split_line(std::string r);

std::vector<std::string> Script(std::string sn);

void seed(unsigned worker_instance);

std::string seed_str();

bool file_exists(const std::string& name);

void find_dump_file(int &suffix);


rapidxml::xml_document<> xml_doc;
rapidxml::xml_node<> * root_node;
rapidxml::xml_node<> * child_node;

std::map<std::string,std::string> configuration;
std::map<std::string,std::string> scripts;

std::vector<std::string> PathwayConfigurations;
double Friction, maxjump_thresh, redo_thresh;
int CoresPerWorker, nRepeats, maxExtraRepeats, globalSeed;
unsigned random_seed;
std::string dump_dir;
bool reseed, seeded, loglammps, postDump, preMin, xml_success,
  spline_path, match_planes,write_dev;

std::map<std::string,std::tuple<double,double,int>> parameters;

private:
  std::mt19937 rng;
};


#endif // PARSER_H
