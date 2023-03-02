#include "Parser.hpp"

std::string lammps_int_to_date(int stamp) {
  std::string month_names[12] =
    {"Jan","Feb","Mar","Apr","May","Jun","Jul","Sep","Oct","Nov","Dec"};

  int year_int = int(stamp / 10000);
  std::string year = std::to_string(year_int);
  int month_int = int(stamp / 100 - year_int*100);
  std::string month = month_names[month_int];
  std::string day = std::to_string(int(stamp % 100));

  return day+" "+month+" "+year;
};

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
	parameters["CoM"] = "1";
  parameters["Friction"] = "0.1";
  parameters["StartCoordinate"] = "0.0";
  parameters["StopCoordinate"] = "1.0";
  parameters["LogLammps"] = "0";
  parameters["MaxJump"] = "0.1";
  parameters["ReSampleThresh"] = "0.5";
  parameters["maxExtraRepeats"] = "1";
  parameters["postMin"] = "0";
	parameters["workerDump"] = "0";
  parameters["PreMin"] = "1";
  parameters["Rediscretize"] = "1";
  parameters["MatchPlanes"] = "0";
  parameters["RealMEPDist"] = "1";
  parameters["FixPAFIGroup"] = "all";
  parameters["CubicSpline"] = "1";


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
  f_error_thresh = std::stod(parameters["ForceErrorThresh"]);
  maxExtraRepeats = std::stoi(parameters["maxExtraRepeats"]);
  postMin = bool(std::stoi(parameters["postMin"]));
	workerDump = bool(std::stoi(parameters["workerDump"]));
  preMin = bool(std::stoi(parameters["PreMin"]));
  rediscretize = bool(std::stoi(parameters["Rediscretize"]));
  use_custom_positions = bool(std::stoi(parameters["UseCustomPositions"]));
  real_coord = int(std::stoi(parameters["RealMEPDist"]));
  cubic_spline = bool(std::stoi(parameters["CubicSpline"]));
  
  std::stringstream ss(parameters["CustomPositions"]);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> tokens(begin, end);
  for (auto &s: tokens) custom_positions.push_back(std::stod(s));
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
  parameters["DumpFolder"] = "./dumps";
  parameters["OverDamped"] = "1";
	parameters["CoM"] = "1";
  parameters["Friction"] = "0.1";
  parameters["StartCoordinate"] = "0.0";
  parameters["StopCoordinate"] = "1.0";
  parameters["LogLammps"] = "0";
  parameters["MaxJump"] = "0.1";
  parameters["ReSampleThresh"] = "0.5";
  parameters["maxExtraRepeats"] = "1";
  parameters["ForceErrorThresh"] = "1.0";
  parameters["UseCustomPositions"] = "0.0";
  parameters["CustomPositions"] = "0.0 1.0";
  std::stringstream ss(parameters["CustomPositions"]);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> tokens(begin, end);
  for (auto &s: tokens) custom_positions.push_back(std::stod(s));
  //parameters["postMin"] = "1";
  //parameters["PreMin"] = "1";
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

void Parser::seed(unsigned _random_seed) {
  random_seed = _random_seed;
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

std::string Parser::seed_str(bool reseed) {
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
  str+="          (c) TD Swinburne and M-C Marinica 2020\n\n";

	str+="\nScripts:\n\n";
	for(auto s: scripts) str+=s.first+" : "+s.second+"\n";

	str+="\nParameters:\n\n";
  for(auto s: parameters) if(s.first!="PathwayConfigurations") str+=s.first+" : "+s.second+"\n";
	str+="\n\n";

  return str;
};

std::vector<double> Parser::sample_r(std::vector<double> pathway_r) {
  // determine knot positions 
  // pathway_r : from NEB
  // All in convention set by RealMEPDist
  std::vector<double> r_vec;
  // Hierarchy: use_custom_positions, rediscretize
  if(use_custom_positions) {
    // overwrite everything with custom positions
    for(auto r: custom_positions) {
      r_vec.push_back(r);
    }
  } else if(rediscretize && nPlanes>1) {
    // Regular interpolation
    double dr = (stopr-startr)/(double)(nPlanes-1);
    for (double r = startr; r <= stopr+0.5*dr; r += dr )
      r_vec.push_back(r);
  } else {
    // Original interpolation
    for(auto r: pathway_r) {
      if(r>=startr-0.02 && r<=stopr+0.02) r_vec.push_back(r);
    }
  }
  return r_vec;
};
  