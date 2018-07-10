#include <stdio.h>
#include "mpi.h"
#include "MumpsSolver.h"

#include <string>

#define USE_COMM_WORLD -987654


class DistributedMatrixEx1
{
public:
  DistributedMatrixEx1(int num_blocks, int block_size, MPI_Comm comm) 
    : comm_(comm), loc_i(NULL), loc_j(NULL), loc_A(NULL), 
      num_blocks_(num_blocks), block_size_(block_size)
  {
    int half_num_blocks = (int) num_blocks/2;
    diagblock_size_ = half_num_blocks*block_size_; // this is the block on the diagonal in the middle of the matrix; it is a diagonal

    n_  = num_blocks*block_size_;
    n_ += diagblock_size_; 
    n_ += block_size*(half_num_blocks-1);

    nnz_  = num_blocks*(block_size_+1)*block_size_/2; //only lower + diagonal entries
    nnz_ += half_num_blocks*block_size_*block_size_;
    nnz_ += diagblock_size_; // this is a diagonal
    //nnz_  += diagblock_size_; // this is the diag block under the above diagonal
    nnz_ += block_size_*(half_num_blocks-1)*4;

    loc_nnz_=0; 

    MPI_Comm_size(comm_, &nranks_);
    MPI_Comm_rank(comm_, &my_rank_);
  }
  virtual ~DistributedMatrixEx1();

  /* distributes the number of nnz equally accross ranks */
  void distributedAssemble_scheme1();
  /* distributes the number of rows across ranks */
  //void distributedAssemble_scheme2();

  /* convert to Fortran indexes */
  void convertToFortran() 
  {
    for (int i=0; i<loc_nnz_; i++) loc_i[i] = loc_i[i]+1;
    for (int i=0; i<loc_nnz_; i++) loc_j[i] = loc_j[i]+1;
  }

  /* outputs the matrix (i,j,value) in a text file: works with one process only */
  void toTextFile() const;
protected:
  MPI_Comm comm_;
  int nranks_;
  int my_rank_;

  int num_blocks_;
  int block_size_;
  int diagblock_size_;
public:
  long long n_;
  long long nnz_, loc_nnz_;
  int* loc_i, *loc_j;
  double* loc_A;


};


int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  DistributedMatrixEx1 dmat(800, 15, MPI_COMM_WORLD);
  dmat.distributedAssemble_scheme1();

  dmat.toTextFile();

  MumpsSolver* mumps = new MumpsSolver(dmat.n_, MPI_COMM_WORLD, MPI_COMM_WORLD);
  dmat.convertToFortran();
  mumps->setLocalEntries(dmat.nnz_, dmat.loc_nnz_, dmat.loc_i, dmat.loc_j, dmat.loc_A);
  mumps->matrixChanged();

  //mumps->saveOrderingPermutation();

  double vec[dmat.n_];
  for(int i=0; i<dmat.n_; i++)
    vec[i]=1.0;
  mumps->solve(vec);

  printf("vec %g %g\n", vec[0], vec[1]);
  delete mumps;

  MPI_Finalize();

  return 0;
}

DistributedMatrixEx1::~DistributedMatrixEx1()
{
  delete[] loc_i; 
  delete[] loc_j;
  delete[] loc_A;
};

void DistributedMatrixEx1::toTextFile() const
{
  if(nranks_!=1) {
    printf("Saving the matrix to a text file works only with one MPI process only\n");
    return;
  }
  std::string file = "matrix.txt";
  FILE *f = fopen(file.c_str(), "w");
  if(NULL!=f) {
    fprintf(f, "%d %d 0\n", n_, nnz_);
    for(long long itnz=0; itnz<nnz_; itnz++) {
      fprintf(f, "%d %d %20.16f\n", loc_i[itnz], loc_j[itnz], loc_A[itnz]);
    }
    fclose(f);
  } else {
    printf("Could not open file %s for writing\n", file.c_str());
  }
}

void DistributedMatrixEx1::distributedAssemble_scheme1()
{
  if(my_rank_==0) printf("Dealing with %d nnz on %d procs\n", nnz_, nranks_);
  //nnz per rank
  long long nnzperrank = nnz_ / nranks_;

  if(my_rank_==0) printf("approx nnz %d per rank\n", nnzperrank);

  long long remainder = nnz_ - nranks_*nnzperrank;
  printf("remainder %d\n", remainder);
  long long loc_nz_start = my_rank_*nnzperrank + (my_rank_<remainder?my_rank_:remainder);
  int next_rank = my_rank_+1;
  long long loc_nz_end = next_rank*nnzperrank + (next_rank<remainder?next_rank:remainder)-1;

  loc_nnz_ = loc_nz_end-loc_nz_start+1;
  printf("proc %d starts at %d-th nnz and ends at %d-th nnz for a total of %d nnz\n", my_rank_, loc_nz_start, loc_nz_end, loc_nnz_);

  loc_i = new int[loc_nnz_];
  loc_j = new int[loc_nnz_];
  loc_A = new double[loc_nnz_];

  int half_num_blocks = (int) num_blocks_/2;

  //build the matrix
  //go over columns.
  long long it_nz=0;
  for(int nb=0; nb<num_blocks_; nb++) {
    
    for(int j=nb*block_size_; j<(nb+1)*block_size_; j++) {
      //diag block
      for(int i=j;i<(nb+1)*block_size_; i++) {
	int nzidx=it_nz-loc_nz_start;
	//is on this rank ?
	if(it_nz>=loc_nz_start && it_nz<=loc_nz_end) {
	  loc_i[nzidx] = i;
	  loc_j[nzidx] = j;
	  loc_A[nzidx] = i==j?2:1.0/block_size_;
	}
	it_nz++;
      }

      //lower subdiag block
      if(nb<half_num_blocks) {
	for(int i=(half_num_blocks+nb)*block_size_; i<(half_num_blocks+nb+1)*block_size_; i++) {
	  if(it_nz>=loc_nz_start && it_nz<=loc_nz_end) {
	    int nzidx=it_nz-loc_nz_start;
	    loc_i[nzidx] = i;
	    loc_j[nzidx] = j;
	    loc_A[nzidx] = 1.0/block_size_;
	  }
	  it_nz++;
	}
      }

      //lower subdiag block - entry 1
      if(nb<half_num_blocks-1) { // for 0 to half_num_blocks-1
	//lower elem per column: (i,j) where i=j + num_blocks_*block_size + half_num_blocks*block_size
	int i=j + num_blocks_*block_size_ + half_num_blocks*block_size_; 
	assert(i<n_); 
	assert(j>=0);
	assert(j<(half_num_blocks-1)*block_size_);
	if(it_nz>=loc_nz_start && it_nz<=loc_nz_end) {
	  int nzidx=it_nz-loc_nz_start;
	  loc_i[nzidx] = i;
	  loc_j[nzidx] = j;
	  loc_A[nzidx] = 1.0/block_size_;
	}
	it_nz++;	
      } 
      //lower subdiag block - entry 2
      if(nb>0 && nb<half_num_blocks) {
	//upper elem per column: (i,j) where i=j + num_blocks_*block_size + (half_num_blocks-1)*block_size
	int i=j + num_blocks_*block_size_ + (half_num_blocks-1)*block_size_;
	assert( i < n_ );
	assert( i >= (half_num_blocks+num_blocks_)*block_size_ );
	assert( j < n_ ); 
	assert( j >= block_size_ );
	if(it_nz>=loc_nz_start && it_nz<=loc_nz_end) {
	  int nzidx=it_nz-loc_nz_start;
	  loc_i[nzidx] = i;
	  loc_j[nzidx] = j;
	  loc_A[nzidx] = 1.0/block_size_;
	}
	it_nz++;
      }
 
    }// end of for j
  } // end of for nb

  // these are the columns for the block on the diagonal in the middle of the matrix; the block is a diagonal and has another diagonal beneath
  for(int j=num_blocks_*block_size_; j<num_blocks_*block_size_+diagblock_size_; j++) {
    if(it_nz>=loc_nz_start && it_nz<=loc_nz_end) {
      int nzidx=it_nz-loc_nz_start;
      loc_i[nzidx] = j;
      loc_j[nzidx] = j;
      loc_A[nzidx] = 1e-6;
    }
    it_nz++;

    int i = j+diagblock_size_;
    if(i>=n_) continue;
    if(it_nz>=loc_nz_start && it_nz<=loc_nz_end) {
      int nzidx=it_nz-loc_nz_start;
      loc_i[nzidx] = i;
      loc_j[nzidx] = j;
      loc_A[nzidx] = 1.0/block_size_;
    }
    it_nz++;
  }

  //these are columns corresponding to the lower diagonal which is size (half_num_blocks-1)*block_size
  for(int j=num_blocks_*block_size_+diagblock_size_; j<n_; j++) {
    assert(j<num_blocks_*block_size_+diagblock_size_ + (half_num_blocks-1)*block_size_);
    if(it_nz>=loc_nz_start && it_nz<=loc_nz_end) {
      int nzidx=it_nz-loc_nz_start;
      loc_i[nzidx] = j;
      loc_j[nzidx] = j;
      loc_A[nzidx] = -1.0;
    }
    it_nz++;
  }
  assert(it_nz==nnz_);

  printf("Created a matrix of size %d with nnz %d (local nnz %d)\n", n_, nnz_, loc_nnz_);

};
