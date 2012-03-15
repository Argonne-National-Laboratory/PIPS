#ifndef SPARSEVECTORHPP
#define SPARSEVECTORHPP

#include "CoinIndexedVector.hpp"
#include "denseVector.hpp"

class sparseVector {
public:
	sparseVector() : len(0) {}
	sparseVector(int len);
	sparseVector(const sparseVector& v2);
	sparseVector(const denseVector& v2);

	~sparseVector() {  }


	bool allocated() const { return (len != 0); }
	int length() const { return len;}
	void deallocate() { len = 0; v.empty(); }
	void allocate(int n) { assert(!allocated()); len = n;  v.reserve(n); }
	// zero out
	void clear();

	//double& operator[](int idx) { assert(idx < len); return v[idx]; }
	const double& operator[](int idx) const { /*assert(idx < len);*/ return v[idx]; }
	
	void divideBy(double pivot);
	void multiplyBy(double x);
	void negate();
	//void etaize(int idx);

	double dotWithSelf() const;
	long double dotWithSelfLongDouble() const;

	void print() const;
	void printDense() const;
	double density() const { return static_cast<double>(v.getNumElements())/len; }


	// just leave this accessable, less wrapping work	
	CoinIndexedVector v;
private:
	int len;

};



#endif
