#include "scenarioReduction.hpp"
#include "scenarioReductionUtilities.hpp"

#include <algorithm>
#include <limits>

using namespace std;

// Fast forward selection, Algorithm 2 in
// "Scenario Reduction and Scenario Tree Construction for Power Management Problems"
// Growe-Kuska, Heitsch, Romisch 2003
// Input: distance matrix, scenario probabilities, total number of scenarios, number of scenarios to pick
// Output: boolean vector, true if scenario is kept, false if scenario is reduced
static void forwardSelection(double const*const* distances, double const *probabilities, bool* keepingScenario, int nScen, int nReduced) {
	
	fill(keepingScenario, keepingScenario+nScen, false);

	vector<int> added;
	added.reserve(nReduced);

	// need to store distances from previous step 
	vector<double> oldC(nScen*nScen), newC(nScen*nScen);
	for (int i = 0; i < nScen; i++)
		for (int j = 0; j < nScen; j++) oldC[i*nScen+j] = distances[i][j];

	for (int i = 0; i < nReduced; i++) {
		int minScen = -1;
		double minVal = numeric_limits<double>::max();
		for (int u = 0; u < nScen; u++) {
			if (keepingScenario[u]) continue;
			double myVal = 0.0;
			for (int k = 0; k < nScen; k++) {
				if (keepingScenario[u] || k == u) continue;
				double val;
				if (i > 0) {
					val = min(oldC[k*nScen+u],oldC[k*nScen+added.back()]);
				} else {
					val = oldC[k*nScen+u];
				}
				newC[k*nScen+u] = val;
				myVal += probabilities[k]*val;
			}
			if (myVal < minVal) {
				minVal = myVal;
				minScen = u;
			}
		}
		newC.swap(oldC);
		assert(minScen >= 0);
		assert(!keepingScenario[minScen]);
		keepingScenario[minScen] = true;
		added.push_back(minScen);
	}

	/*
	oldC.clear(); newC.clear();
	vector<vector<distanceScenarioPair> > sortedDistances;
	sortDistances(distances, nScen, sortedDistances);

	double obj = calculateObjective(sortedDistances, probabilities, keepingScenario);
	printf("Objective for fast forward selection: %f\n",obj);
	*/

}



vector<int> fastForwardSelection(stochasticInput &input, int nScenariosWanted) {
	vector<int> out;
	if (nScenariosWanted == 0) return out;

	double **distances;
	generateDistances(distances,input);
	int nScen = input.nScenarios();
	
	bool *keepingScenario = new bool[nScen];
	vector<double> probs(nScen);
	for (int i = 0; i < nScen; i++) probs[i] = input.scenarioProbability(i);


	forwardSelection(distances, &probs[0], keepingScenario, nScen, nScenariosWanted);

	freeDistances(distances,input);

	for (int i = 0; i < nScen; i++) {
		if (keepingScenario[i]) out.push_back(i);
	}
	
	delete [] keepingScenario;
	
	assert(out.size() == static_cast<unsigned>(nScenariosWanted));
	
	return out;
}
