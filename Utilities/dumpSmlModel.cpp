#include "sml.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>

using namespace std;

/*
Dump raw problem data to files.
Format assumes that only bounds change per scenario.

First-stage file format:
<number of first-stage variables> <number of first-stage constraints>
<number of scenarios>
Only bounds vary
<number of second-stage variables> <number of second-stage constraints>
<for each first-stage variable, a row with: <name> <lower bound> <upper bound> <objective coefficient>>
<for each first-stage constraint, a row with: <name> <lower bound> <upper bound>>
A matrix
<number nonzero in A matrix>
<A matrix column starts>
<A matrix row indices>
<A matrix entries>
W matrix
<number nonzero in W matrix>
<W matrix column starts>
<W matrix row indices>
<W matrix entries>
T matrix
<number nonzero in T matrix>
<T matrix column starts>
<T matrix row indices>
<T matrix entries>

Second-stage file format:
<number of second-stage variables> <number of second-stage constraints>
Only Bounds Vary
<for each second-stage variable, a row with: <name> <lower bound> <upper bound> <objective coefficient>>
<for each second-stage constraint, a row with: <name> <lower bound> <upper bound>>


See the code that reads/writes these files if this isn't clear.
*/


void dumpSmlModelOnlyBoundsVary(const string &modelfilename, const string &datafilename, const string &outbasename) {
	
	ExpandedModelInterface *emroot = sml_generate(modelfilename, datafilename, false);
	ExpandedModelInterface *em = emroot->children.at(0);
	
	int nScenarios = em->children.size();
	int nFirstStageVar = em->getNLocalVars(), nFirstStageCons = em->getNLocalCons();
	int nSecondStageVar = em->children[0]->getNLocalVars(), nSecondStageCons = em->children[0]->getNLocalCons();

	stringstream f0name;
	f0name << outbasename << "0";
	ofstream f0(f0name.str().c_str());
	f0 << nFirstStageVar << " " << nFirstStageCons << endl;
	f0 << nScenarios << endl;
	f0 << "Only Bounds Vary\n";
	f0 << nSecondStageVar << " " << nSecondStageCons << endl;

	vector<double> lb(nFirstStageVar), ub(nFirstStageVar), obj(nFirstStageVar);
	em->getObjGradient(&obj[0]);
	em->getColLowBounds(&lb[0]);
	em->getColUpBounds(&ub[0]);
	
	const list<string> &colnames = em->getLocalVarNames();
	list<string>::const_iterator it = colnames.begin();
	f0.precision(15);
	for (int i = 0; i < nFirstStageVar; ++i, ++it) {
		f0 << *it << " " << lb[i] << " " << ub[i] << " " << obj[i] << endl; 
	}
	assert(it == colnames.end());

	lb.resize(nFirstStageCons); ub.resize(nFirstStageCons);
	em->getRowBounds(&lb[0],&ub[0]);
	const list<string> &rownames = em->getLocalConNames();
	it = rownames.begin();
	for (int i = 0; i < nFirstStageCons; ++i, ++it) {
		f0 << *it << " " << lb[i] << " " << ub[i] << endl; 
	}
	assert(it == rownames.end());

	f0 << "A matrix\n";

	{
		int nnz = em->getNzJacobianOfIntersection(em);
		vector<int> colbeg(nFirstStageVar+1);
		vector<int> collen(nFirstStageVar);
		vector<int> rowidx(nnz);
		vector<double> elts(nnz);

		em->getJacobianOfIntersection(em, &colbeg[0], &collen[0], &rowidx[0], &elts[0]);
		f0 << nnz << endl;
		for (int i = 0; i <= nFirstStageVar; i++) {
			if (i != nFirstStageVar) assert(colbeg[i] + collen[i] == colbeg[i+1]);
			if (i != 0) f0 << " ";
			f0 << colbeg[i];
		}
		f0 << endl;
		for (int i = 0; i < nnz; i++) {
			assert(rowidx[i] < nFirstStageCons);
			if (i != 0) f0 << " ";
			f0 << rowidx[i];
		}
		if (nnz) f0 << endl;
		for (int i = 0; i < nnz; i++) {
			if (i != 0) f0 << " ";
			f0 << scientific << elts[i];
		}
		if (nnz) f0 << endl;
	}

	f0 << "W matrix\n";
	{
		int nnz = em->children[0]->getNzJacobianOfIntersection(em->children[0]);
		vector<int> colbeg(nSecondStageVar+1);
		vector<int> collen(nSecondStageVar);
		vector<int> rowidx(nnz);
		vector<double> elts(nnz);

		em->children[0]->getJacobianOfIntersection(em->children[0], &colbeg[0], &collen[0], &rowidx[0], &elts[0]);
		f0 << nnz << endl;
		for (int i = 0; i <= nSecondStageVar; i++) {
			if (i != nSecondStageVar) assert(colbeg[i] + collen[i] == colbeg[i+1]);
			if (i != 0) f0 << " ";
			f0 << colbeg[i];
		}
		f0 << endl;
		for (int i = 0; i < nnz; i++) {
			assert(rowidx[i] < nSecondStageCons);
			if (i != 0) f0 << " ";
			f0 << rowidx[i];
		}
		f0 << endl;
		for (int i = 0; i < nnz; i++) {
			if (i != 0) f0 << " ";
			f0 << scientific << elts[i];
		}
		f0 << endl;
	}

	f0 << "T matrix\n";
	{
		int nnz = em->children[0]->getNzJacobianOfIntersection(em);
		vector<int> colbeg(nFirstStageVar+1);
		vector<int> collen(nFirstStageVar);
		vector<int> rowidx(nnz);
		vector<double> elts(nnz);

		em->children[0]->getJacobianOfIntersection(em, &colbeg[0], &collen[0], &rowidx[0], &elts[0]);
		f0 << nnz << endl;
		for (int i = 0; i <= nFirstStageVar; i++) {
			if (i != nFirstStageVar) assert(colbeg[i] + collen[i] == colbeg[i+1]);
			if (i != 0) f0 << " ";
			f0 << colbeg[i];
		}
		f0 << endl;
		for (int i = 0; i < nnz; i++) {
			assert(rowidx[i] < nSecondStageCons);
			if (i != 0) f0 << " ";
			f0 << rowidx[i];
		}
		f0 << endl;
		for (int i = 0; i < nnz; i++) {
			if (i != 0) f0 << " ";
			f0 << scientific << elts[i];
		}
		f0 << endl;
	}

	f0.close();

	for (int k = 0; k < nScenarios; k++) {
		stringstream f1name;
		f1name << outbasename << k+1;
		ofstream f1(f1name.str().c_str());

		f1 << nSecondStageVar << " " << nSecondStageCons << endl;
		f1 << "Only Bounds Vary\n";

		lb.resize(nSecondStageVar); ub.resize(nSecondStageVar); obj.resize(nSecondStageVar);
		em->children[k]->getObjGradient(&obj[0]);
		em->children[k]->getColLowBounds(&lb[0]);
		em->children[k]->getColUpBounds(&ub[0]);
		
		const list<string> &colnames = em->children[k]->getLocalVarNames();
		it = colnames.begin();
		f1.precision(15);
		for (int i = 0; i < nSecondStageVar; ++i, ++it) {
			f1 << *it << " " << lb[i] << " " << ub[i] << " " << obj[i] << endl; 
		}
		assert(it == colnames.end());

		lb.resize(nSecondStageCons); ub.resize(nSecondStageCons);
		em->children[k]->getRowBounds(&lb[0],&ub[0]);

		const list<string> &rownames = em->children[k]->getLocalConNames();
		it = rownames.begin();
		f1.precision(15);
		for (int i = 0; i < nSecondStageCons; ++i, ++it) {
			f1 << *it << " " << lb[i] << " " << ub[i] << endl; 
		}
		assert(it == rownames.end());

		f1.close();
		// this will mess with internal SML destructor, but we need to free the memory
		delete em->children[k];
		em->children[k] = 0;


	}

	delete emroot;

}


int main(int argc, char **argv) {
	if (argc != 4) {
		printf("usage: %s modelfile datafile outbasename\n",argv[0]);
		return 0;
	}
	dumpSmlModelOnlyBoundsVary(argv[1],argv[2],argv[3]);

	return 0;
}
