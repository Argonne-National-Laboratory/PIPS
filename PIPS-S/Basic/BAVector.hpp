#ifndef BAVECTOR_HPP
#define BAVECTOR_HPP

#include "BA.hpp"
#include "sparseVector.hpp"
#include "denseVector.hpp"
#include "memContainer.hpp"
#include <cmath>

//#include <iostream>

const int MAX_UPDATES = 500;

// would be nice to use "enum class" from C++11 here
enum BAVectorType { DualVector, PrimalVector, BasicVector };

// Basic is for holding RHS vectors for FTRAN and BTRAN
// see note at BALinearAlgebra1::ftran


// T1 and T2 indicate vector types to use for first and second stage vectors, respectively
// these should be either sparseVector or denseVector
// with some changes, maybe you could get a tree strucure (of static depth) by nesting this class with itself
template<typename T1, typename T2> class BAVector {
public:

  BAVector() : dims(0), nScen(0) {
    
  }
private:
  BAVector(const BAVector&) { }
  BAVector<T1,T2>& operator=(const BAVector<T1,T2>&) { }
public:
	BAVector(const BADimensions &dims, const BAContext &ctx, BAVectorType t) {
		allocate(dims,ctx,t);
	}

	~BAVector() {
		for (int i = 0; i < nScen; i++) {
			if (vec2[i]) delete vec2[i];
		}
	}

	void allocate(const BADimensions &dims, const BAContext &ctx, BAVectorType t) {
	  assert(this->dims == 0);
	  vecType = t;
	  nScen = dims.numScenarios();
	  vec2.resize(nScen);
	  this->dims = &dims;
	  this->ctx = &ctx;
	  localScen.reserve(10); 
	  localScen.push_back(-1);
	  if (vec1.allocated()) vec1.deallocate();
	  if (t == PrimalVector) {
	    vec1.allocate(dims.numFirstStageVars());
	  } else if (t == DualVector) {
	    vec1.allocate(dims.numFirstStageCons());
	  } else { // Basic
	    vec1.allocate(dims.numFirstStageVars()+MAX_UPDATES); // basic vars >= cons
	  }
	  for (int i = 0; i < nScen; i++) {
	    if (vec2[i]) delete vec2[i];
	    if (ctx.assignedScenario(i)) {
	      vec2[i] = new T2();
	      if (t == PrimalVector) {
		vec2[i]->allocate(dims.numSecondStageVars(i));
	      } else if (t == DualVector) {
		vec2[i]->allocate(dims.numSecondStageCons(i));
	      } else { // Basic
		vec2[i]->allocate(dims.numSecondStageCons(i)+MAX_UPDATES); // cons >= basic vars
	      }
	      localScen.push_back(i);
	    } else {
	      vec2[i] = 0;
	    }
	  }
	  
	  
	}
  
  bool allocated() const { return (dims != 0); }
  
  void divideBy(double d) {
    assert(vecType == BasicVector); 
    vec1.divideBy(d);
    for (int i = 0; i < nScen; i++) {
      vec2[i]->divideBy(d);
    }
  }
  
	void multiplyBy(double d) {
		assert(vecType == BasicVector);
		vec1.multiplyBy(d);
		for (int i = 0; i < nScen; i++) {
			vec2[i]->multiplyBy(d);
		}
	}

	// this could use some rethinking
	// if we're basic and call negate, do we actually
	// want to negate all vectors or only local vectors?
	void negate() {
		vec1.negate();
		for (int i = 0; i < nScen; i++) {
			if (vec2[i]) vec2[i]->negate();
		}
	}

	void clear() {
		vec1.clear();
		for (int i = 0; i < nScen; i++) {
			if (vec2[i]) vec2[i]->clear();
		}
	}


	double dotWithSelf() const {
		double d = 0;
		for (unsigned i = 1; i < localScen.size(); i++) {
			d += vec2[localScen[i]]->dotWithSelf();
		}
		d = ctx->reduce(d);
		d += vec1.dotWithSelf();
		return d;
	}

	T1& getFirstStageVec() { return vec1; }
	const T1& getFirstStageVec() const { return vec1; }

	T2& getSecondStageVec(int i) { assert(i >= 0 && vec2[i]); return *vec2[i]; }
	const T2& getSecondStageVec(int i) const { assert(i >= 0 && vec2[i]); return *vec2[i]; }

	bool hasScenario(int i) const { if (i == -1) return true; else return (vec2[i] != 0); }
	const std::vector<int>& localScenarios() const { return localScen; }

	BAVectorType vectorType() const { return vecType; }

  
protected:
	const BADimensions* dims;
	int nScen;
	T1 vec1;
	// NULL means not assigned this scenario
	std::vector<T2*> vec2;
	BAVectorType vecType;
	std::vector<int> localScen;
	const BAContext *ctx;


};

// special accessors when first and second stage are the same type
template<typename T> class sameBAVector : public BAVector<T,T> {
public:
  sameBAVector() : BAVector<T,T>() {}
  sameBAVector(const BADimensions &dims, const BAContext &ctx, BAVectorType t)
    : BAVector<T,T>(dims,ctx,t) {}
private:
private:
  sameBAVector(const sameBAVector&) { }
  sameBAVector<T>& operator=(const sameBAVector<T>&) { }
public:
	// so we can iterate through first and second stages uniformly
	T& getVec(int idx) {
		if (idx == -1) return this->vec1;
		else { assert(this->vec2[idx]); return *this->vec2[idx]; }
	}

	const T& getVec(int idx) const {
		if (idx == -1) return this->vec1;
		else { assert(this->vec2[idx]); return *this->vec2[idx]; }
	}


};

class denseBAVector : public sameBAVector<denseVector> {
public:
  denseBAVector() : sameBAVector<denseVector>() {}
  //denseBAVector(const BADimensions &dims, const BAContext &ctx, BAVectorType t) : BAVector<denseVector,denseVector>(dims,ctx,t) {}
  void copyFrom(const denseBAVector &v) {
    vec1.copyFrom(v.getFirstStageVec());
    for (unsigned i = 1; i < localScen.size(); i++) {
      int idx = localScen[i];
      vec2[idx]->copyFrom(v.getSecondStageVec(idx));
    }
  }
  
 void copyBeginning(const denseBAVector &v) {
    vec1.copyBeginning(v.getFirstStageVec());
    for (unsigned i = 1; i < localScen.size(); i++) {
      int idx = localScen[i];
      vec2[idx]->copyBeginning(v.getSecondStageVec(idx));
    }
  }

  denseBAVector(const denseBAVector& v)
    : sameBAVector<denseVector>()
  {
    allocate(*v.dims,*v.ctx,v.vecType);
    copyFrom(v);
  }		 

	
  denseBAVector& operator=(const denseBAVector& v) {
    //these are required by "allocate" to work 
    //also the dealocation of the local vectors is done in allocate
    dims=0; nScen=0; localScen.clear();
    allocate(*v.dims,*v.ctx,v.vecType);
    copyFrom(v);
  }

	double dotWith(const denseBAVector &v) {
		double dot = 0;
		for (unsigned i = 1; i < localScen.size(); i++) {
			int idx = localScen[i];
			dot += vec2[idx]->dotWith(v.getSecondStageVec(idx));
		}
		dot = ctx->reduce(dot);
		dot += vec1.dotWith(v.getFirstStageVec());
		return dot;
	}
	
	inline double& operator[](BAIndex i) {
		if (i.scen == -1) {
			return this->vec1[i.idx];
		} else {
			return (*this->vec2[i.scen])[i.idx];
		}
	}
	inline double operator[](BAIndex i) const {
		if (i.scen == -1) {
			return this->vec1[i.idx];
		} else {
			return (*this->vec2[i.scen])[i.idx];
		}
	}


};


class sparseBAVector : public sameBAVector<sparseVector> {
public:
	sparseBAVector()  {}

private:
sparseBAVector(const sparseBAVector&) {}
  sparseBAVector& operator=(const sparseBAVector&) {}
public:
	inline double operator[](BAIndex i) const {
		if (i.scen == -1) {
			return this->vec1[i.idx];
		} else {
			return (*this->vec2[i.scen])[i.idx];
		}
	}
	

	void checkClean() {
		assert(vecType == BasicVector);
		for (int scen = -1; scen < nScen; scen++) {
			getVec(scen).v.checkClean();
		}

	}

protected:
	// operator[] gets confused with const/non-const. try this for now
	inline double& operator()(BAIndex i) {
		if (i.scen == -1) {
			return this->vec1.v.denseVector()[i.idx];
		} else {
			return this->vec2[i.scen]->v.denseVector()[i.idx];
		}
	}
};


// distributed vector for any type
template<typename T> class BAFlagVector {
public:
	BAFlagVector(const BADimensions &dims, const BAContext &ctx, BAVectorType t) : dims(0),vec1(0) {
		allocate(dims, ctx, t);
	}

	BAFlagVector() : dims(0), nScen(0), vec1(0) {}

  BAFlagVector(const BAFlagVector& v) : dims(0),vec1(0), nScen(0) {
    this->allocate(*v.dims, *v.ctx, v.vecType);
    this->copyFrom(v);
  }

  BAFlagVector<T>& operator=(const BAFlagVector<T>& v) {
    if(this==&v) return *this;
    this->deallocate();
    this->allocate(*v.dims, *v.ctx, v.vecType);
    this->copyFrom(v);
    return *this;
  }  

	~BAFlagVector() {
		deallocate();	
	}

	void allocate(const BADimensions &dims, const BAContext &ctx, BAVectorType t) {
		assert(this->dims == 0);
		this->dims = &dims;
		this->ctx = &ctx;
		vecType = t;
		nScen = dims.numScenarios();
		vec2.resize(nScen);
		if (t == PrimalVector) {
			vec1 = new denseFlagVector<T>(dims.numFirstStageVars());
		} else {
			vec1 = new denseFlagVector<T>(dims.numFirstStageCons());
		}
		for (int i = 0; i < nScen; i++) {
			if (ctx.assignedScenario(i)) {
				if (t == PrimalVector) {
					vec2[i] = new denseFlagVector<T>(dims.numSecondStageVars(i));
				} else {
					vec2[i] = new denseFlagVector<T>(dims.numSecondStageCons(i));
				}
			} else {
				vec2[i] = 0;
			}
		}
	}

	void deallocate() {
		if (vec1) delete vec1;
		for (unsigned i = 0; i < vec2.size(); i++) {
			if (vec2[i]) delete vec2[i];
		}
		dims = 0;
	}


	denseFlagVector<T>& getFirstStageVec() { return *vec1; }
	const denseFlagVector<T>& getFirstStageVec() const { return *vec1; }

	denseFlagVector<T>& getSecondStageVec(int i) { assert(vec2[i]); return *vec2[i]; }
	const denseFlagVector<T>& getSecondStageVec(int i) const { assert(vec2[i]); return *vec2[i]; }

	// shouldn't use this for looping, only for isolated accesses
	inline T& operator[](BAIndex i) {
		if (i.scen == -1) {
			return (*vec1)[i.idx];
		} else {
			return (*vec2[i.scen])[i.idx];
		}
	}

	inline const T& operator[](BAIndex i) const {
		if (i.scen == -1) {
			return (*vec1)[i.idx];
		} else {
			return (*vec2[i.scen])[i.idx];
		}
	}

	denseFlagVector<T>& getVec(int i) {
		if (i == -1) { return *vec1; }
		else { return *vec2[i]; }
	}

	const denseFlagVector<T>& getVec(int i) const {
		if (i == -1) { return *vec1; }
		else return *vec2[i];
	}
	BAVectorType vectorType() const { return vecType; }


	bool hasScenario(int i) const { if (i == -1) return true; else return (vec2[i] != 0); }

	void copyFrom(const BAFlagVector<T> &r) {
		vec1->copyFrom(r.getFirstStageVec());
		for (int scen = 0; scen < nScen; scen++) {
			if (vec2[scen]) vec2[scen]->copyFrom(r.getSecondStageVec(scen));
		}
	
	}

protected:
	const BADimensions *dims;
	int nScen;
	denseFlagVector<T> *vec1;
	std::vector<denseFlagVector<T>*> vec2;
	BAVectorType vecType;
	const BAContext *ctx;
};

// a bit different from BAVector.
// used to store list of basic variables and vector of primal infeasibilities
// T must be a container type, like std::vector or CoinIndexedVector
template<typename T> class BAContainer {
public:
	

	BAContainer() : dims(0) {}

	BAContainer(const BADimensions &dims, const BAContext &ctx, BAVectorType t) : dims(0) {
		allocate(dims,ctx,t);
	}
  BAContainer(const BAContainer& v) : dims(0) {
    allocate(*v.dims, *v.ctx, v.vecType);
    vec1=v.vec1;
    vec2=v.vec2;
  }

  BAContainer<T>& operator=(const BAContainer<T>& v) {
    if(this==&v) return *this;
    deallocate();
    allocate(*v.dims, *v.ctx, v.vecType);

    vec1=v.vec1;
    vec2=v.vec2;

    return *this;
  }
  


	void allocate(const BADimensions &dims, const BAContext &ctx, BAVectorType t) {
		//assert(this->dims == 0);
		vecType = t;
		nScen = dims.numScenarios();
		vec2.resize(nScen);
		this->dims = &dims;
		this->ctx = &ctx;
		localScen.clear();
		localScen.push_back(-1);
		if (t == PrimalVector) {
			vec1.reserve(dims.numFirstStageVars());
		} else if (t == DualVector) {
			vec1.reserve(dims.numFirstStageCons());
		} else { // Basic
		  vec1.reserve(dims.numFirstStageVars()+MAX_UPDATES); // basic vars >= cons
		}
		for (int i = 0; i < nScen; i++) {
			if (ctx.assignedScenario(i)) {
				localScen.push_back(i);
				if (t == PrimalVector) {
					vec2[i].reserve(dims.numSecondStageVars(i));
				} else if (t == DualVector) {
					vec2[i].reserve(dims.numSecondStageCons(i));
				} else { // Basic
					vec2[i].reserve(dims.numSecondStageCons(i)+MAX_UPDATES);
				}
			} 
		}
	}

  void deallocate() {
    clear(); 
    dims = 0;
  }
	void clear() {
		vec1.clear();
		for (unsigned i = 1; i < localScen.size(); i++) {
			vec2[localScen[i]].clear();
		}
	}

	T& getFirstStageVec() { return vec1; }
	const T& getFirstStageVec() const { return vec1; }

	T& getSecondStageVec(int i) { assert(i >= 0); return vec2[i]; }
	const T& getSecondStageVec(int i) const { assert(i >= 0); return vec2[i]; }

	const std::vector<int>& localScenarios() const { return localScen; }

	T& getVec(int i) {
		if (i == -1) return vec1;
		else { return vec2[i]; }
	}

	const T& getVec(int i) const {
		if (i == -1) return vec1;
		else { return vec2[i]; }
	}


protected:
	const BADimensions* dims;
	int nScen;
	T vec1;
	std::vector<T> vec2;
	BAVectorType vecType;
	std::vector<int> localScen;
  const BAContext* ctx;
};

// second iteration of etaBAVector
// storage for local parts of packed eta vectors
// routines to apply the eta-matrices aren't here
class packedBAVector {
public:
	packedBAVector() {}

	void allocate(const BADimensions &dims, const BAContext &ctx) {
		vec.resize(ctx.localScenarios().size());
	}

	// takes a local index (index for a localScenarios() array)
	// see usage in BAPFIPAR.cpp
	packedVector& getPackedVec(int localidx) {
		return vec.at(localidx);
	}

	const packedVector& getPackedVec(int localidx) const {
		return vec.at(localidx);
	}
	
	// release allocated memory
	void clear() {
		for (unsigned i = 0; i < vec.size(); i++) {
			vec[i].resizeAndDestroy(0);
		}
	}

protected:
	std::vector<packedVector> vec;

};



#endif

