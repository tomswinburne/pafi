#include "Parser.hpp"

// Loads config file, lammps scripts etc

Parser::Parser(std::string file, bool test) {
  xml_success = false;

  // Default Values simulation independent...
  configuration["CoresPerWorker"]="1";

  configuration["SampleSteps"] = "1000";
  configuration["ThermSteps"] = "1000";
  configuration["ThermWindow"] = "500";
  configuration["nRepeats"] = "1";
  configuration["DumpFolder"] = "./dumps";
  configuration["OverDamped"] = "1";
  configuration["Friction"] = "0.1";
  configuration["LogLammps"] = "0";
  configuration["MaxJump"] = "0.1";
  configuration["ReSampleThresh"] = "0.5";
  configuration["maxExtraRepeats"] = "1";
  configuration["PostDump"] = "0";
  configuration["PreMin"] = "1";
  configuration["SplinePath"] = "1";
  configuration["MatchPlanes"] = "0";
  configuration["GlobalSeed"] = "137";
  configuration["FreshSeed"] = "1";
  // these could be overwritten during TI, or set to the lambda=1 value
  configuration["LinearThermalExpansionX"] = "0.0";
  configuration["LinearThermalExpansionY"] = "0.0";
  configuration["LinearThermalExpansionZ"] = "0.0";
  configuration["QuadraticThermalExpansionX"] = "0.0";
  configuration["QuadraticThermalExpansionY"] = "0.0";
  configuration["QuadraticThermalExpansionZ"] = "0.0";


  seeded = false;
  // Read the xml file into a vector
	std::ifstream xmlfile(file);
	std::vector<char> buffer((std::istreambuf_iterator<char>(xmlfile)),
    std::istreambuf_iterator<char>());
	buffer.push_back('\0');

  // Parse the buffer using the xml file parsing library into xml_doc
	xml_doc.parse<0>(&buffer[0]);

	root_node = xml_doc.first_node();

  bool found_setup = false, found_scripts = false, found_scan = false;

  while (root_node) {

    if(rtws(root_node->name())=="Setup") {
      found_setup = true;
      child_node = root_node->first_node();
      while (child_node) {
        if(rtws(child_node->name())=="PathwayConfigurations") {
          PathwayConfigurations = split_lines(child_node->value());
        } else {
          configuration[rtws(child_node->name())] = rtws(child_node->value());
        }
        child_node = child_node->next_sibling();
      }
    } else if(rtws(root_node->name())=="Parameters") {
      found_scan = true;
      child_node = root_node->first_node();
      while (child_node) {
        std::vector<std::string> ss = split_line(child_node->value());
        if(ss.size()!=3)
          std::cout<<"Error in "<<child_node->value()<<" XML declaration"<<std::endl;
        parameters[rtws(child_node->name())] =
          std::make_tuple(std::stod(ss[0]),std::stod(ss[1]),std::stoi(ss[2]));

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

  if(!found_setup) {
    std::cout<<"XML file incomplete! Provide Configuration!"<<std::endl;
    return;
  }
  if(!found_scan) {
    std::cout<<"XML file incomplete! Provide Parameters!"<<std::endl;
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
  CoresPerWorker = std::stoi(configuration["CoresPerWorker"]);
	nRepeats = std::stoi(configuration["nRepeats"]);

  int therm_window = std::max(1,std::stoi(configuration["ThermSteps"])/2);
  configuration["ThermWindow"] = std::to_string(therm_window);

  dump_dir = configuration["DumpFolder"];
	Friction = std::stod(configuration["Friction"]);
	loglammps = bool(std::stoi(configuration["LogLammps"]));
  maxjump_thresh = std::stod(configuration["MaxJump"]);
  redo_thresh = std::stod(configuration["ReSampleThresh"]);
  maxExtraRepeats = std::stoi(configuration["maxExtraRepeats"]);
  postDump = bool(std::stoi(configuration["PostDump"]));
  preMin = bool(std::stoi(configuration["PreMin"]));
  spline_path = bool(std::stoi(configuration["SplinePath"]));
  match_planes = !bool(std::stoi(configuration["Rediscretize"]));
  globalSeed = std::stoi(configuration["GlobalSeed"]);
  reseed = bool(std::stoi(configuration["FreshSeed"]));
  write_dev = bool(std::stoi(configuration["WriteDev"]));
};

void Parser::insert_params(std::string &s, Holder &params) {
  std::string s_r;
  for(auto param : params) {
    std::regex e ("%"+param.first+"%");
    s_r = std::regex_replace(s,e,std::to_string(param.second));
    s = s_r;
  }
};

void Parser::overwrite_xml(int nProcs) {
  // Default Values
  configuration["CoresPerWorker"]=std::to_string(nProcs);
  configuration["SampleSteps"] = "1";
  configuration["ThermSteps"] = "1";
  configuration["ThermWindow"] = "1";
  configuration["nRepeats"] = "1";
  configuration["DumpFolder"] = "./dumps";
  configuration["OverDamped"] = "1";
  configuration["Friction"] = "0.1";
  configuration["MaxJump"] = "1.0";
  configuration["ReSampleThresh"] = "0.5";
  configuration["maxExtraRepeats"] = "0";
  configuration["PostDump"] = "1";
  configuration["PreMin"] = "1";
  configuration["WriteDev"] = "0";
  parameters["Temperature"] = std::make_tuple(0.,0.,1);
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

std::vector<std::string> Parser::split_line(std::string r) {
  std::istringstream iss(r);
  std::vector<std::string> f{std::istream_iterator<std::string>{iss},
    std::istream_iterator<std::string>{}};
  return f;
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
	//configuration["Temperature"] = std::to_string(T);
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

	str+="\nConfiguration:\n\n";
  for(auto s: configuration) if(s.first!="PathwayConfigurations") str+=s.first+" : "+s.second+"\n";
	str+="\n\n";

  str+="\nParameter Scans :\n\n";
  std::string tab = "\t";
  for(auto s= std::prev(parameters.end()); ; s=std::prev(s)) {
    str += tab;
    str += s->first+" : "+std::to_string(std::get<0>(s->second));
    str += " -> "+std::to_string(std::get<1>(s->second))+" in ";
    str += std::to_string(std::get<2>(s->second))+" steps\n\n";
    tab += "\t";
    if(s==parameters.begin()) break;
  }
  str+="\n";
  return str;
};

// Check for dump files. Ugly implementation for portability across filesystems
void Parser::find_dump_file(int &suffix) {
  std::ofstream raw;
  std::string params_file;
  FILE *file;
  for (suffix=0; suffix < 100; suffix++) {
    params_file = dump_dir+"/params_"+std::to_string(suffix);
    file = fopen(params_file.c_str(), "r"); // read only - don't overwrite !
    if (file!=NULL) {
        fclose(file);
        continue;
    }
    raw.open(params_file.c_str(),std::ofstream::out);
    if(raw.is_open()) {
      raw<<welcome_message();
      raw.close();
    }
    break;
    if(suffix==100) suffix=-1;
  }
};
