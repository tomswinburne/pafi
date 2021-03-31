#include "Parser.hpp"

// Loads config file, lammps scripts etc

Parser::Parser(std::string file, bool test) {
  xml_success = false;

  // Default Values
  parameters["CoresPerWorker"]="1";
  parameters["LowTemperature"] = "100";
  parameters["HighTemperature"] = "1000";
  parameters["TemperatureSteps"] = "10";
  parameters["LinearThermalExpansionX"] = "0.0";
  parameters["LinearThermalExpansionY"] = "0.0";
  parameters["LinearThermalExpansionZ"] = "0.0";
  parameters["QuadraticThermalExpansionX"] = "0.0";
  parameters["QuadraticThermalExpansionY"] = "0.0";
  parameters["QuadraticThermalExpansionZ"] = "0.0";
  parameters["SampleSteps"] = "1000";
  parameters["ThermSteps"] = "1000";
  parameters["ThermWindow"] = "500";
  parameters["nRepeats"] = "1";
  parameters["nPlanes"] = "10";
  parameters["DumpFolder"] = "./dumps";
  parameters["OverDamped"] = "1";
  parameters["Friction"] = "0.1";
  parameters["StartCoordinate"] = "0.0";
  parameters["StopCoordinate"] = "1.0";
  parameters["LogLammps"] = "0";
  parameters["MaxJump"] = "0.1";
  parameters["ReSampleThresh"] = "0.5";
  parameters["maxExtraRepeats"] = "1";
  parameters["PostDump"] = "0";
  parameters["PreMin"] = "1";
  parameters["SplinePath"] = "1";
  parameters["MatchPlanes"] = "0";
  parameters["GlobalSeed"] = "137";
  parameters["FreshSeed"] = "1";


  seeded = false;
  // Read the xml file into a vector
	std::ifstream xmlfile(file);
	std::vector<char> buffer((std::istreambuf_iterator<char>(xmlfile)), std::istreambuf_iterator<char>());
	buffer.push_back('\0');

  // Parse the buffer using the xml file parsing library into xml_doc
	xml_doc.parse<0>(&buffer[0]);

	root_node = xml_doc.first_node();

  bool found_pafi = false, found_scripts = false;

  while (root_node) {

    if(rtws(root_node->name())=="PAFI") {
      found_pafi = true;
      child_node = root_node->first_node();
      while (child_node) {
        if(rtws(child_node->name())=="PathwayConfigurations") {
          PathwayConfigurations = split_lines(child_node->value());
        } else {
          parameters[rtws(child_node->name())] = rtws(child_node->value());
        }
        child_node = child_node->next_sibling();
      }
    } else if(rtws(root_node->name())=="Scripts") {
      found_scripts = true;
      child_node = root_node->first_node();
      while (child_node) {
        scripts[rtws(child_node->name())] = child_node->value();
        child_node = child_node->next_sibling();
      }
    }
    root_node = root_node->next_sibling();
  }
  if(!found_pafi) {
    std::cout<<"XML file incomplete! Provide PAFI parameters!"<<std::endl;
    return;
  }

  if(!found_scripts) {
    std::cout<<"XML file incomplete! Provide MD scripts!"<<std::endl;
    return;
  }

  xml_success = true;

  if(!test) set_parameters();

};

void Parser::set_parameters() {
  // Now we can convert to type
  CoresPerWorker = std::stoi(parameters["CoresPerWorker"]);
	nPlanes = std::stoi(parameters["nPlanes"]);
	nRepeats = std::stoi(parameters["nRepeats"]);
  int therm_window = std::max(1,std::stoi(parameters["ThermSteps"])/2);
  parameters["ThermWindow"] = std::to_string(therm_window);
  dump_dir = parameters["DumpFolder"];
	lowT = std::stod(parameters["LowTemperature"]);
	highT = std::stod(parameters["HighTemperature"]);
	Friction = std::stod(parameters["Friction"]);
	TSteps = std::stoi(parameters["TemperatureSteps"]);
  startr = std::stod(parameters["StartCoordinate"]);
  stopr = std::stod(parameters["StopCoordinate"]);
  loglammps = bool(std::stoi(parameters["LogLammps"]));
  maxjump_thresh = std::stod(parameters["MaxJump"]);
  redo_thresh = std::stod(parameters["ReSampleThresh"]);
  maxExtraRepeats = std::stoi(parameters["maxExtraRepeats"]);
  postDump = bool(std::stoi(parameters["PostDump"]));
  preMin = bool(std::stoi(parameters["PreMin"]));
  spline_path = bool(std::stoi(parameters["SplinePath"]));
  match_planes = !bool(std::stoi(parameters["Rediscretize"]));
  globalSeed = std::stoi(parameters["GlobalSeed"]);
  reseed = bool(std::stoi(parameters["FreshSeed"]));
};

void Parser::overwrite_xml(int nProcs) {
  // Default Values
  parameters["CoresPerWorker"]=std::to_string(nProcs);
  parameters["LowTemperature"] = "0";
  parameters["HighTemperature"] = "0";
  parameters["TemperatureSteps"] = "1";
  parameters["SampleSteps"] = "1";
  parameters["ThermSteps"] = "1";
  parameters["ThermWindow"] = "1";
  parameters["nRepeats"] = "1";
  //parameters["nPlanes"] = "10";
  parameters["DumpFolder"] = "./dumps";
  parameters["OverDamped"] = "1";
  parameters["Friction"] = "0.1";
  parameters["StartCoordinate"] = "0.0";
  parameters["StopCoordinate"] = "1.0";
  parameters["LogLammps"] = "0";
  parameters["MaxJump"] = "0.1";
  parameters["ReSampleThresh"] = "0.5";
  parameters["maxExtraRepeats"] = "1";
  parameters["PostDump"] = "1";
  parameters["PreMin"] = "1";

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

void Parser::seed(unsigned worker_instance) {
  random_seed = (worker_instance+1)*globalSeed;
  rng.seed(random_seed);
	seeded=true;
};

std::vector<std::string> Parser::split_lines(std::string r) {
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
    random_seed=0;
		seeded=true;
	}
  unsigned r = random_seed;
  // give exactly the same seed each time unless reseed=True
  if(reseed) {
    std::uniform_int_distribution<unsigned> d(1,1000000);
    r=d(rng);
  }
  return std::to_string(r);
};

std::vector<std::string> Parser::Script(std::string sn) {
	//parameters["Temperature"] = std::to_string(T);
	return split_lines(scripts[sn]);
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
  str+="          (c) TD Swinburne and M-C Marinica 2021\n\n";

	str+="\nScripts:\n\n";
	for(auto s: scripts) str+=s.first+" : "+s.second+"\n";

	str+="\nParameters:\n\n";
  for(auto s: parameters) if(s.first!="PathwayConfigurations") str+=s.first+" : "+s.second+"\n";
	str+="\n\n";

  return str;
};

/*
Check for dump files. Ugly implementation for portability across filesystems
*/
bool Parser::file_exists(const std::string& name) {
    FILE *file = fopen(name.c_str(), "r");
    if (file!=NULL) {
        fclose(file);
        return true;
    } else {
        return false;
    }
};

void Parser::find_dump_file(std::ofstream &raw, int &suffix){
  std::string params_file;
  for (suffix=0; suffix < 100; suffix++) {
    params_file = dump_dir+"/params_"+std::to_string(suffix);
    if(!file_exists(params_file)) {
      raw.open(params_file.c_str(),std::ofstream::out);
      if(raw.is_open()) {
        raw<<welcome_message();
        raw.close();
        break;
      }
    }
  }
  if(suffix==100) suffix=-1;
};
