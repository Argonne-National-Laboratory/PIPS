#ifndef SCENREDUTILS_HPP
#define SCENREDUTILS_HPP

#include "stochasticInput.hpp"

struct distanceScenarioPair{ 
	double d; int scen; 
	bool operator<(const distanceScenarioPair &s) const { return (d < s.d); } 
};
// sort distances per scenario, saving the original indices
// outer vector i contains sorted vector of distances between other scenarios and scenario i
void sortDistances(double const*const* distances, int nScen, std::vector<std::vector<distanceScenarioPair> > &out);

// calculates the objective to the optimal reduction problem
double calculateObjective(const std::vector<std::vector<distanceScenarioPair> > &sortedDistances,  
		double const *probabilities, bool const *keepingScenario);

// generate distance matrix
// caller is responsible for freeing, using convenience function below
void generateDistances(double **&distances, stochasticInput&);
template <typename T> void freeDistances(T **d, stochasticInput&in) {
	int nscen = in.nScenarios();
	for (int i = 0; i < nscen; i++) delete [] d[i];
	delete [] d;
}


#endif
