#ifndef DENSEVECTORHPP
#define DENSEVECTORHPP

#include <algorithm>
#include <cassert>

class denseVector {
public:
	denseVector() : d(0), own(true), len(0) {}
	denseVector(int len) : own(true), len(len) { d = new double[len]; }
	//denseVector(int len, double *buf) : d(buf), own(false), len(len) {}
	denseVector(const denseVector& v) : own(true), len(v.len) { 
		d = new double[len];
		std::copy(v.d,v.d+len,d);
	}

	~denseVector() { if (d && own) delete [] d; }


	bool allocated() const { return (d != 0); }
	int length() const { assert(d); return len;}
	void deallocate() { assert(d); assert(own); delete [] d; d = 0; len = 0;}
	void allocate(int n) { assert(d == 0); d = new double[n]; len = n; }
	// zero out
	void clear() { assert(d); std::fill(d,d+len,0.0); }

	inline double& operator[](int idx) { /*assert(d); assert(idx < len);*/ return d[idx]; }
	inline const double& operator[](int idx) const { /*assert(d); assert(idx < len);*/ return d[idx]; }
private:
  denseVector& operator=(const denseVector&) {}
public:
	void divideBy(double pivot);
	void multiplyBy(double x);
	void negate();
	void etaize(int idx);

	double dotWith(const denseVector &eta) const;
	double dotWithSelf() const;
	void add(const denseVector &eta);

	void print() const;
	void copyFrom(const double *b,int n);
	void copyFrom(const denseVector &);
	void copyBeginning(const double *b,int n);
	void copyBeginning(const denseVector &v, int n);
	void copyToPosition(const double *b, int start, int n);

	void copyToPosition(const denseVector &v, int start, int n);
void copyBeginning(const denseVector &v);
	
	void swap(denseVector &v) {
		assert(v.own == own);
		int l = v.len; double *d2 = v.d;
		v.len = len; v.d = d;
		len = l; d = d2;
	}

	  double* getPointer() { return d; }

private:
	double *d;
	const bool own; // own buffer?
	int len;

};

template <class T> class denseFlagVector {
	
public:
	denseFlagVector(int len) : len(len) { d = new T[len]; }
private:
  denseFlagVector(const denseFlagVector&) {}
  denseFlagVector& operator=(const denseFlagVector&) {}
public:
	virtual ~denseFlagVector() { delete [] d; }

	virtual int length() const { return len; }
        virtual bool allocated() const { return (d != 0); }
  
	inline T& operator[](int idx) { /*assert(idx < len);*/ return d[idx]; }
	inline const T& operator[](int idx) const { /*assert(idx < len);*/ return d[idx]; }
	void copyFrom(const denseFlagVector<T> &v) {
		assert(v.len == len);
		for (int i = 0; i < len; i++) d[i] = v.d[i];
	}

	void copyBeginning(const denseFlagVector<T> &v, int n){
	  assert(allocated() && len >= v.length());
		copy(v.d,v.d+v.length(),d);
	}

	void copyToPosition(const T *b, int start, int n){
		assert(allocated() && len >= start+n);
		copy(b,b+n,&d[start]);
	}

private:
	T *d;
	const int len;

};

#endif
