#ifndef BAPFIPAR_HPP
#define BAPFIPAR_HPP

#include "BAData.hpp"

/* Manager for parallel PFI updates
  doesn't include FTRAN/BTRAN calls as this is in template wrapper BAPFIParWrapper

  Store packed (in-order) column-ETAs, in addition to dense
  vector for each eta containing the entries in pivot rows. This dense vector is used
  for BTRAN-PFI and "mini" FTRAN-PFI
*/


class BAPFIPar {
public:
	BAPFIPar(const BAData &d, int max_updates = MAX_UPDATES);
	~BAPFIPar();



	// rhs vector is result of FTRAN
	void ftranPFI(sparseBAVector &rhs);

	void btranPFI(sparseBAVector &rhs) { assert(npfi == 0); }

	// will pull leaving index from each eta and add to pivotEntries
	// must be called before btranPFISimplex
	void setLeaving(BAIndex leaving);

	// assumes that rhs is zero except for a single unit entry
	// corresponding to variable that has been selected to leave
	// CANNOT handle arbitrary rhs
	void btranPFISimplex(sparseBAVector &rhs);

	// add a new eta vector. ftranVec is result of FTRAN
	void newEta(sparseBAVector &ftranVec, BAIndex in, BAIndex out);
	
	// clear ETAs
	void clear();

	int nUpdates() const { return npfi; }

private:
	
	// do mini-ftran to determine pivot values by a dense triangular solve
	// rhs should contain "leaving" entries from ftran RHS
	// on exit, will contain sequence of pivot values
	void miniFtranDense(CoinIndexedVector &rhs) const;
	void miniFtranSparsish(CoinIndexedVector &rhs) const;

	const BAData &data;
	std::vector<BAIndex> enter;
	std::vector<BAIndex> leave;
	std::vector<packedBAVector> etas;
	/* column-oriented matrix of pivot entries
	   in (i,j), entry in leave[i] from eta[j]
	   this is actually a strictly lower-triangular matrix
	   with the convention that pivot entry in etas are stored in *entering* position
	   also entries are NEGATED from traditional form so that dense triangular solve can be performed
	*/
	double *pivotEntries, *pivotRow;
	std::vector<std::vector<int> > pivotEntriesIdx; // column-oriented list of nonzeros in pivotEntries
	CoinIndexedVector pivotRhs; // work buffer for mini ftran and btran
	int npfi, nScen;
	const int maxUpdates;
	BAIndex leaving;

	// localLeaving[p][i] == k means that process p has scenario leave[k].scen
	// used for collecting rhs for mini ftran
	std::vector<std::vector<int> > localLeaving;
	// firstStageLeaving[i] == k means leave[k].scen == -1 (so every process has it)
	std::vector<int> firstStageLeaving;
	unsigned maxLocalLeaving; // for Allgather

	memContainer<double> ftranRecv;

	std::vector<int> scenToLocalIdx;

	// not safe to copy
	BAPFIPar(const BAPFIPar&);
	BAPFIPar& operator=(const BAPFIPar&);

};





#endif
