#include "Parser.hpp"

// Loads config file, lammps scripts etc

Parser::Parser(std::string file) {
  seeded = false;
  // Read the xml file into a vector
	std::ifstream xmlfile(file);
	std::vector<char> buffer((std::istreambuf_iterator<char>(xmlfile)), std::istreambuf_iterator<char>());
	buffer.push_back('\0');

  // Parse the buffer using the xml file parsing library into xml_doc
	xml_doc.parse<0>(&buffer[0]);
	root_node = xml_doc.first_node("PAFI");

  parameters["CoresPerWorker"]=\
    rtws(root_node->first_node("CoresPerWorker")->value());
  parameters["LowTemperature"] = \
    rtws(root_node->first_node("LowTemperature")->value());
  parameters["HighTemperature"] = \
    rtws(root_node->first_node("HighTemperature")->value());
  parameters["TemperatureSteps"] = \
    rtws(root_node->first_node("TemperatureSteps")->value());
  parameters["LinearThermalExpansion"] = \
    rtws(root_node->first_node("LinearThermalExpansion")->value());
  parameters["QuadraticThermalExpansion"] = \
    rtws(root_node->first_node("QuadraticThermalExpansion")->value());
  parameters["SampleSteps"] = \
    rtws(root_node->first_node("SampleSteps")->value());
  parameters["ThermSteps"] = \
    rtws(root_node->first_node("ThermSteps")->value());
  parameters["ThermWindow"] = \
    rtws(root_node->first_node("ThermWindow")->value());
	parameters["nRepeats"] = \
    rtws(root_node->first_node("nRepeats")->value());
  parameters["nPlanes"] = \
    rtws(root_node->first_node("nPlanes")->value());
  parameters["DumpFolder"] = \
    rtws(root_node->first_node("DumpFolder")->value());
  parameters["OverDamped"] = \
    rtws(root_node->first_node("OverDamped")->value());
  parameters["Friction"] = \
    rtws(root_node->first_node("Friction")->value());
  parameters["StartCoordinate"] = \
    rtws(root_node->first_node("StartCoordinate")->value());
  parameters["StopCoordinate"] = \
    rtws(root_node->first_node("StopCoordinate")->value());
  parameters["LogLammps"] = \
    rtws(root_node->first_node("LogLammps")->value());

	// Now we can convert to type
	KnotList = Parse(root_node->first_node("KnotList")->value());

  CoresPerWorker = std::stoi(parameters["CoresPerWorker"]);
	nPlanes = std::stoi(parameters["nPlanes"]);
	nRepeats = std::stoi(parameters["nRepeats"]);
	//dump_dir = std::to_string(parameters["DumpFolder"]);
  dump_dir = parameters["DumpFolder"];
	lowT = std::stod(parameters["LowTemperature"]);
	highT = std::stod(parameters["HighTemperature"]);
	Friction = std::stod(parameters["Friction"]);
	TSteps = std::stoi(parameters["TemperatureSteps"]);
  startr = std::stod(parameters["StartCoordinate"]);
  stopr = std::stod(parameters["StopCoordinate"]);
  loglammps = bool(std::stoi(parameters["LogLammps"]));

  rapidxml::xml_node<> * scrs_n = root_node->first_node("Scripts");

  for (rapidxml::xml_node<> * scr_n = scrs_n->first_node("Script"); scr_n; scr_n = scr_n->next_sibling()) {
    std::string key = scr_n->first_attribute("name")->value();
    scripts[key] = scr_n->value();
  }

};

// remove leading and trailing whitespaces - avoids RapidXML parsing bug
std::string Parser::rtws(std::string s) {
	int bws=0,ews=0;
	for(auto c=s.begin();c!=s.end();c++,bws++)
		if(!std::isspace(*c)) break;
	for(auto c=std::prev(s.end());c!=s.begin();c--,ews++)
		if(!std::isspace(*c)) break;
	return s.substr(bws,s.length()-ews-bws);
};

void Parser::seed(int random_seed) {
  rng.seed(random_seed);
	seeded=true;
};

std::vector<std::string> Parser::Parse(std::string r) {
  std::vector< std::string > sc;
  std::istringstream input;
  input.str(r);
  for (std::string line; std::getline(input, line); ) {
    std::string nline = rtws(line);
    if(!nline.empty()) sc.push_back(nline);
  }
	return sc;
};

std::string Parser::seed_str() {
  if(!seeded) {
		std::cout<<"NOT SEEDED!!!\n";
		rng.seed(0);
		seeded=true;
	}
  std::uniform_int_distribution<unsigned> d(1,1000000);
  unsigned r=d(rng);
  return std::to_string(r);
};

std::vector<std::string> Parser::Script(std::string sn) {
	//parameters["Temperature"] = std::to_string(T);
	return Parse(scripts[sn]);
};

// Pointless :)

std::string Parser::welcome_message(){
  std::string str;
	str="\n";
  str+="       _______      _______      _______     _________\n";
  str+="      (  ____ )    (  ___  )    (  ____ \\    \\__   __/\n";
  str+="      | (    )|    | (   ) |    | (    \\/       ) (\n";
  str+="      | (____)|    | (___) |    | (__           | |\n";
  str+="      |  _____)    |  ___  |    |  __)          | |\n";
  str+="      | (          | (   ) |    | (             | |\n";
  str+="      | )          | )   ( |    | )          ___) (___\n";
  str+="      |/           |/     \\|    |/           \\_______/\n";
  str+="      Projected    Average      Force        Integrator\n";
  str+="          (c) TD Swinburne and M-C Marinica 2018\n\n";

	str+="\nScripts:\n\n";
	for(auto s: scripts) str+=s.first+" : "+s.second+"\n";

	str+="\nParameters:\n\n";
  for(auto s: parameters) if(s.first!="KnotList") str+=s.first+" : "+s.second+"\n";
	str+="\n\n";

  return str;
};
