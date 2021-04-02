
#define BOLTZ 8.617e-5
#define BAR2EVA3 6.25e-7

typedef std::map<std::string,double> Holder;


template<typename A,typename B>
using PairHolder = std::map<std::string,std::pair<A,B>>;
/*
template<typename A>
using VectorHolder = std::map<std::string,std::vector<A>> ;
//std::map<Holder,VectorHolder<double>> all_results;
*/
