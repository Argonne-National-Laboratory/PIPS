#ifndef BAPFIPARWRAP_HPP
#define BAPFIPARWRAP_HPP

// wrap linear algebra class LA with PFI-providing class PFI

template <typename LA, typename PFI> class BAPFIParWrapper {
public:
	BAPFIParWrapper(const BAData& d) : la(d), pfi(d) {}
	~BAPFIParWrapper() {}

	void reinvert(const BAFlagVector<variableState> &s) {
		la.reinvert(s);
		pfi.clear();
	}
	void ftran(sparseBAVector &rhs) {
		la.ftran(rhs);
		pfi.ftranPFI(rhs);
	}
	void btran(sparseBAVector &rhs) {
		pfi.btranPFI(rhs); // may only work if just after reinversion
		la.btran(rhs);
	}

	// this is btran as used in plain simplex
	// rhs must be a unit vector with 1-entry in position leaving
	// cannot call multiple times between updates (as currently implemented)
	void btranSimplex(sparseBAVector &rhs, BAIndex leaving) {
		pfi.setLeaving(leaving);
		pfi.btranPFISimplex(rhs);
		la.btran(rhs);
	}

	// ftranVec is ftran applied to entering column
	int replaceColumnNoF(BAIndex enter, BAIndex leave, sparseBAVector &ftranVec) {
		pfi.newEta(ftranVec,enter,leave);
		return 0; // no stability checks here
	}

	// only useful for testing purposes
	int replaceColumn(BAIndex enter, BAIndex leave, const sparseBAVector &vec) {
		ftran(vec);
		return replaceColumnNoF(enter,leave,vec);
	}

	int nUpdates() const { return pfi.nUpdates(); }

protected:
	LA la;
	PFI pfi;
};


#endif

