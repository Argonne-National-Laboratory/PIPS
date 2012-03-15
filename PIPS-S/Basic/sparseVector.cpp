#include "sparseVector.hpp"

#include <iostream>

using namespace std;

sparseVector::sparseVector(int len) : v(len), len(len) {}


sparseVector::sparseVector(const sparseVector& v2) : v(v2.v), len(v2.len) { }

sparseVector::sparseVector(const denseVector& v2) : len(v2.length()) {
	v.reserve(len);
	int nonZero = 0;
	int *idx = v.getIndices();
	double *elts = v.denseVector();
	for (int i = 0; i < len; i++) {
		if (v2[i]) {
			idx[nonZero++] = i;
			elts[i] = v2[i];
		}
	}
	v.setNumElements(nonZero);
}

void sparseVector::divideBy(double pivot) {
	v /= pivot;
}

void sparseVector::multiplyBy(double x) {
	v *= x;
}

void sparseVector::negate() {
	double *elts = v.denseVector();
	const int *idx = v.getIndices();
	int nelts = v.getNumElements();
	for (int i = 0; i < nelts; i++) {
		int r = idx[i];
		elts[r] = -elts[r];
	}
}

void sparseVector::clear() {
	v.clear();
}

/*
void sparseVector::etaize(int idx) {
	double pivot = v[idx];
	pivot = 1/pivot;
	multiplyBy(pivot);
	negate();
	v[idx] = pivot;
}*/



void sparseVector::print() const {
	int nelt = v.getNumElements();
	if (nelt == 0) {
		cout << "No elements\n";
	} else {
		const int *idx = v.getIndices();
		const double *elts = v.denseVector();
		for (int i = 0; i < nelt; i++) {
			cout << idx[i] << ": " << elts[idx[i]] << endl;	
		}
	}
}

void sparseVector::printDense() const {
	for (int i = 0; i < len; i++) {
		cout << v[i] << " ";
	} cout << endl;

}


double sparseVector::dotWithSelf() const {
	double dot = 0;
	int nelt = v.getNumElements();
	const int *idxs = v.getIndices();
	const double *elts = v.denseVector();

	for (int i = 0; i < nelt; i++) {
		int idx = idxs[i];
		dot += elts[idx]*elts[idx];
	}
	return dot;
}

long double sparseVector::dotWithSelfLongDouble() const {
	long double dot = 0;
	int nelt = v.getNumElements();
	const int *idxs = v.getIndices();
	const double *elts = v.denseVector();
	for (int i = 0; i < nelt; i++) {
		int idx = idxs[i];
		long double val = elts[idx];
		dot += val*val;
	}
	return dot;
}
