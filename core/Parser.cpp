#include "Parser.hpp"

// Loads config file, lammps scripts etc

Parser::Parser(std::string file) {

  seeded = false;

  // Parse the XML into the property tree.
  boost::property_tree::read_xml(file,tree,\
  boost::property_tree::xml_parser::no_comments);

  // | boost::property_tree::xml_parser::trim_whitespace); // RapidXML can bug
  // => have to send everything through home rolled rtws(std::string s)

  parameters["CoresPerWorker"]=\
    rtws(tree.get<std::string>("PAFI.CoresPerWorker","1"));
  parameters["LowTemperature"] = \
    rtws(tree.get<std::string>("PAFI.LowTemperature","0."));
  parameters["HighTemperature"] = \
    rtws(tree.get<std::string>("PAFI.HighTemperature","100."));
  parameters["TemperatureSteps"] = \
    rtws(tree.get<std::string>("PAFI.TemperatureSteps","2"));
  parameters["LinearThermalExpansion"] = \
    rtws(tree.get<std::string>("PAFI.LinearThermalExpansion","0.0"));
  parameters["QuadraticThermalExpansion"] = \
    rtws(tree.get<std::string>("PAFI.QuadraticThermalExpansion","0.0"));
  parameters["SampleSteps"] = \
    rtws(tree.get<std::string>("PAFI.SampleSteps","100"));
  parameters["SampleWindow"] = \
    rtws(tree.get<std::string>("PAFI.SampleWindow","100"));
  parameters["ThermSteps"] = \
    rtws(tree.get<std::string>("PAFI.ThermSteps","100"));
    parameters["ThermWindow"] = \
      rtws(tree.get<std::string>("PAFI.ThermWindow","100"));
  parameters["nPlanes"] = \
    rtws(tree.get<std::string>("PAFI.nPlanes","100"));
  parameters["DumpFolder"] = \
    rtws(tree.get<std::string>("PAFI.DumpFolder","dumps"));
  parameters["OverDamped"] = \
    rtws(tree.get<std::string>("PAFI.OverDamped","1"));
  parameters["Friction"] = \
    rtws(tree.get<std::string>("PAFI.Friction","0.05"));
  parameters["StartCoordinate"] = \
    rtws(tree.get<std::string>("PAFI.StartCoordinate","0.0"));
  parameters["StopCoordinate"] = \
    rtws(tree.get<std::string>("PAFI.StopCoordinate","1.0"));

	// Now we can convert to type
	KnotList = Parse(tree.get<std::string>("PAFI.KnotList"),false);

  CoresPerWorker = boost::lexical_cast<int>(parameters["CoresPerWorker"]);
	nPlanes = boost::lexical_cast<int>(parameters["nPlanes"]);
	dump_dir = boost::lexical_cast<std::string>(parameters["DumpFolder"]);
	lowT = boost::lexical_cast<double>(parameters["LowTemperature"]);
	highT = boost::lexical_cast<double>(parameters["HighTemperature"]);
	Friction = boost::lexical_cast<double>(parameters["Friction"]);
	TSteps = boost::lexical_cast<int>(parameters["TemperatureSteps"]);
  startr = boost::lexical_cast<double>(parameters["StartCoordinate"]);
  stopr = boost::lexical_cast<double>(parameters["StopCoordinate"]);


	BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("PAFI.Scripts")) {
		std::string key =  boost::lexical_cast<std::string>(v.first);
		scripts[key] = tree.get<std::string>("PAFI.Scripts."+v.first);
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
}

std::vector<std::string> Parser::Parse(std::string r, bool needseed) {
	if(!seeded && needseed) {
		std::cout<<"NOT SEEDED!!!\n";
		rng.seed(0);
		seeded=true;
	}

  std::string raw=r;
	//replacements based on parameters
	for(auto it=parameters.begin(); it!=parameters.end(); it++ ) {
		std::string key=it->first;
		key="%"+key+"%";
		boost::trim(key);
		boost::replace_all(raw, key, it->second);
	}
  // random seed
  boost::random::uniform_01<> uniform;
	boost::random::uniform_int_distribution<> d(1,1000000);
  std::string key="%RANDOM%";
	while(not boost::find_first(raw, key).empty()) {
		int r=d(rng);
		std::string s=boost::lexical_cast<std::string>(boost::format("%1%" ) % r );
		boost::replace_first(raw, key, s);
	}
  // Split Lines and remove trailing whitespace
  std::vector< std::string > sc,ssc;
	boost::split( sc, raw, boost::is_any_of("\n"), boost::token_compress_off );
	for(auto s:sc){
		std::string ns = rtws(s);
		if(!ns.empty()) ssc.push_back(ns);
	}
	return ssc;
};

std::vector<std::string> Parser::Script(std::string sn) {
	//parameters["Temperature"] = boost::lexical_cast<std::string>(T);
	return Parse(scripts[sn]);
};

// Pointless :)
void Parser::welcome_message(){
	std::cout<<"\n";
  std::cout<<"       _______      _______      _______     _________\n";
  std::cout<<"      (  ____ )    (  ___  )    (  ____ \\    \\__   __/\n";
  std::cout<<"      | (    )|    | (   ) |    | (    \\/       ) (\n";
  std::cout<<"      | (____)|    | (___) |    | (__           | |\n";
  std::cout<<"      |  _____)    |  ___  |    |  __)          | |\n";
  std::cout<<"      | (          | (   ) |    | (             | |\n";
  std::cout<<"      | )          | )   ( |    | )          ___) (___\n";
  std::cout<<"      |/           |/     \\|    |/           \\_______/\n";
  std::cout<<"      Projected    Average      Force        Integrator\n";
  std::cout<<"          (c) TD Swinburne and M-C Marinica 2018\n\n";

	std::cout<<"\nScripts:\n\n";

	for(auto s: scripts) std::cout<<s.first<<" : "<<s.second<<"\n";

	std::cout<<"\nParameters:\n\n";

	BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("PAFI")) if(v.first != "Scripts" && v.first != "KnotList") {
		std::cout<<"\t"<<v.first<<" : "<<parameters[v.first]<<"\n";
	}
	std::cout<<"\n\n";

};
