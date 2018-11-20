#ifndef PARSER_H
#define PARSER_H

// remove trailing whitespace - avoids RapidXML parsing bug
std::string rtws(std::string s) {
	int bws=0,ews=0;
	for(auto c=s.begin();c!=s.end();c++,bws++)
		if(!std::isspace(*c)) break;
	for(auto c=std::prev(s.end());c!=s.begin();c--,ews++)
		if(!std::isspace(*c)) break;
	return s.substr(bws,s.length()-ews-bws);
}

// Loads config file, lammps scripts etc
class Parser {
public:
Parser(std::string file) {
	boost::random::random_device rd;
	rng.seed(rd());

	// Parse the XML into the property tree.
	boost::property_tree::read_xml(file,tree,\
    boost::property_tree::xml_parser::no_comments);
    // | boost::property_tree::xml_parser::trim_whitespace);



	parameters["CoresPerWorker"]=\
		rtws(tree.get<std::string>("PAFI.CoresPerWorker","1"));

  parameters["LowTemperature"] = \
		rtws(tree.get<std::string>("PAFI.LowTemperature","100."));

	parameters["HighTemperature"] = \
		rtws(tree.get<std::string>("PAFI.HighTemperature","100."));

	parameters["position"] = \
		rtws(tree.get<std::string>("PAFI.Position","-1."));

	parameters["Linear"] = \
		rtws(tree.get<std::string>("PAFI.LinearThermalExpansion","0.0"));

	parameters["Quadratic"] = \
		rtws(tree.get<std::string>("PAFI.QuadraticThermalExpansion","0.0"));

  parameters["SampleSteps"] = \
		rtws(tree.get<std::string>("PAFI.SampleSteps","100"));

	parameters["ThermSteps"] = \
	 rtws(tree.get<std::string>("PAFI.ThermSteps","100"));

	parameters["TWindow"] = \
 	 rtws(tree.get<std::string>("PAFI.ThermWindow","100"));

	parameters["nPlanes"] = rtws(tree.get<std::string>("PAFI.nPlanes","100"));

	KnotList = Parse(tree.get<std::string>("PAFI.KnotList"));

  CoresPerWorker = boost::lexical_cast<int>(parameters["CoresPerWorker"]);
	nPlanes = boost::lexical_cast<int>(parameters["nPlanes"]);

	BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("PAFI.Scripts")) {
		std::string key =  boost::lexical_cast<std::string>(v.first);
		scripts[key] = tree.get<std::string>("PAFI.Scripts."+v.first);
	}

};


std::vector<std::string> Parse(std::string r) {
  boost::random::mt11213b rng;
  boost::random::random_device rd;
	rng.seed(rd());

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

std::vector<std::string> Script(std::string sn) {
	//parameters["Temperature"] = boost::lexical_cast<std::string>(T);
	return Parse(scripts[sn]);
};

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
