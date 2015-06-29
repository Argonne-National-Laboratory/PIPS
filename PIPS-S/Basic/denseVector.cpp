#include "denseVector.hpp"

#include "PIPSLogging.hpp"

#include <iostream>
#include <cstdio>

using namespace std;

// could call BLAS for some of these, but not worth the trouble

void denseVector::divideBy(double pivot) {
	for (int i = 0; i < len; i++) {
		d[i] /= pivot;
	}
}

void denseVector::multiplyBy(double x) {
	for (int i = 0; i < len; i++) {
		d[i] *= x;
	}
}

void denseVector::negate() {
	for (int i = 0; i < len; i++) {
		d[i] = -d[i];
	}
}

void denseVector::etaize(int idx) {
	double pivot = d[idx];
	pivot = 1/pivot;
	multiplyBy(pivot);
	negate();
	d[idx] = pivot;
}

void denseVector::add(const denseVector &v) {

	for (int i = 0; i < len; i++) {
		d[i] += v[i];
	}
}

double denseVector::dotWith(const denseVector &v) const{

	double dot = 0;
	for (int i = 0; i < len; i++) {
		dot += d[i]*v[i];
	}
	return dot;
}


void denseVector::print() const {
	assert(d);
	for (int i = 0; i < len; i++) {
		cout << d[i] << " ";	
	}
	cout << endl;
}

void denseVector::copyFrom(const double *b, int n) {
	if (allocated() && len != n) {
		assert(0);
		deallocate(); 
		PIPS_APP_LOG_SEV(warning) << "wrong size input buffer: "
					  << "source has length " << n
					  << ", dest has length " << len;
		allocate(n);
	} else if (!allocated()) {
		allocate(n);
	}
	copy(b,b+n,d);
}

void denseVector::copyFrom(const denseVector &v) {
	int n = v.length();
	if (allocated() && len != n) {
		assert(0);
		deallocate(); 
		PIPS_APP_LOG_SEV(warning) << "wrong size input buffer: "
					  << "source has length " << n
					  << ", dest has length " << len;
		allocate(n);
	} else if (!allocated()) {
		allocate(n);
	}
	copy(v.d,v.d+n,d);
}

void denseVector::copyBeginning(const double *b, int n) {
	assert(allocated() && len >= n);
	copy(b,b+n,d);
}

double denseVector::dotWithSelf() const {
	double dot = 0;
	for (int i = 0; i < len; i++) {
		dot += d[i]*d[i];
	}
	return dot;
}
