#ifndef MEMCONTAINER_HPP
#define MEMCONTAINER_HPP

// very simple lightweight memory container
// unlike std::vector, doesn't initialize elements, doesn't copy elements if resized

template<typename T> class memContainer {
public:
	memContainer() : capacity(0), mem(0) {}
	~memContainer() { if (mem) delete [] mem; }

	// don't allow copying if initialized
	memContainer(const memContainer& m) : capacity(0), mem(0) { assert(m.mem == 0); }
	memContainer& operator=(const memContainer& m) { assert(m.mem == 0); return *this; }

	void resizeAndDestroy(int cap) {
		if (cap > capacity) {
			if (mem) delete [] mem;
			mem = new T[cap+100];
			capacity = cap+100;
		} else if (cap == 0) {
			if (mem) delete [] mem;
			mem = 0; capacity = 0;
		}
	}

	const T* ptr() const { return mem; }
	T* ptr() { return mem; }
	int getCapacity() const { return capacity; }

protected:
	int capacity;
	T *mem;
};

// very simple packed vector. follow CoinIndexedVector API where appropriate
// meant to serve as persistent storage for ETA vectors (see BAPFIPar)

class packedVector {
public:	
	packedVector() : len(0) {}

	int *getIndices() { return idx.ptr(); } 
	const int *getIndices() const { return idx.ptr(); }

	double *denseVector() { return vals.ptr(); }
	const double *denseVector() const { return vals.ptr(); }

	void resizeAndDestroy(int n) {
		idx.resizeAndDestroy(n);
		vals.resizeAndDestroy(n);
	}
	int getCapacity() const { return idx.getCapacity(); }
	void setNumElements(int n) { assert(n <= getCapacity()); len = n; }
	int getNumElements() const { return len; }



protected:
	memContainer<int> idx;
	memContainer<double> vals;
	int len;
};



#endif
