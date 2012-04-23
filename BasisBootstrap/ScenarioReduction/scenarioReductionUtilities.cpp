#include "scenarioReductionUtilities.hpp"
#include <algorithm>
#include <cmath>

using namespace std;

void sortDistances(double const*const* distances, int nScen, vector<vector<distanceScenarioPair> > &out) {
out.resize(nScen);
	for (int i = 0; i < nScen; i++) {
		out[i].resize(nScen-1);
		int pos = 0;
		for (int j = 0; j < nScen; j++) {
			if (i == j) continue;
			out[i][pos].d = distances[i][j];
			out[i][pos++].scen = j;
		}
		sort(out[i].begin(),out[i].end());
	}
}

double calculateObjective(const vector<vector<distanceScenarioPair> > &sortedDistances,  
		double const *probabilities, bool const *keepingScenario) {
	double obj = 0.0;
	int nScen = sortedDistances.size();
	for (int i = 0; i < nScen; i++) {
		if (keepingScenario[i]) continue;
		int j;
		for (j = 0; !keepingScenario[sortedDistances[i][j].scen]; j++);
		obj += probabilities[i]*sortedDistances[i][j].d;
	}
	return obj;	

}
namespace{
double vectorDiff2(vector<double> const &v1, vector<double> const &v2) {
	double sum = 0;
	for (unsigned i = 0; i < v1.size(); i++) {
		double diff = v1[i]-v2[i];
		sum += diff*diff;
	}
	return sum;

}
double calculateDistance(stochasticInput &input, int s1, int s2) {

	assert(input.scenarioDimensionsEqual());
	assert(input.onlyBoundsVary());
	// skipping second-stage objective
	vector<double> const &l1 = input.getSecondStageColLB(s1), 
		l2 = input.getSecondStageColLB(s2),
		u1 = input.getSecondStageColUB(s1), 
		u2 = input.getSecondStageColUB(s2), 
		bl1 = input.getSecondStageRowLB(s1), 
		bl2 = input.getSecondStageRowLB(s2),
		bu1 = input.getSecondStageRowUB(s1),
		bu2 = input.getSecondStageRowUB(s2);

	double d = vectorDiff2(l1,l2) + vectorDiff2(u1,u2) + 
		vectorDiff2(bl1,bl2) + vectorDiff2(bu1,bu2);
	return sqrt(d);
}
}
void generateDistances(double **&distances, stochasticInput& input) {
	int nscen = input.nScenarios();
	distances = new double*[nscen];
	for (int i = 0; i < nscen; i++) distances[i] = new double[nscen];

	for (int i = 0; i < nscen; i++)
		for (int j = 0; j < nscen; j++) distances[i][j] = calculateDistance(input,i,j);
	
}
