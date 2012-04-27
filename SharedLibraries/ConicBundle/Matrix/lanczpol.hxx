/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/lanczpol.hxx

    This file is part of ConciBundle, a C/C++ library for convex optimization.

    ConicBundle is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ConicBundle is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

***************************************************************************** */



#ifndef CH_MATRIX_CLASSES__LANCZPOL_HXX
#define CH_MATRIX_CLASSES__LANCZPOL_HXX

/**  @file lanczpol.hxx
    @brief Header declaring the abstract classes CH_Matrix_Classes::Lanczosmatrix and CH_Matrix_Classes::Lanczos 
    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg

*/


#ifndef CH_MATRIX_CLASSES__LANCZOS_HXX
#include "lanczos.hxx"
#endif
#include "clock.hxx"

namespace CH_Matrix_Classes {


/**@addtogroup Lanczos_Interface
*/
  //@{ 


  /** @brief A Lanczos method allowing spectral transformation by Chebycheff polynomials and premature termination

  The code is a translation and adaptation of a FORTRAN code most likely
  written by Hua. 

  */

class Lanczpol:public Lanczos,protected Memarrayuser
{
private:
  int     ierr;        ///< error return code 
  Integer maxop;       ///< upper bound on matrix vector multiplications
  Integer maxiter;     ///< upper bound on number of restarts
  Integer maxguessiter;///< upper bound on number of restarts to guess spectral interval
  Integer guessmult;   ///< number of matrix vector multiplications to guess interval
  Integer choicencheb; ///< user's choice for number of block Chebychev iterations (<0 -> automatic determination)
  Integer nchebit;     ///< number of block Chebychev iterations within one iteration
  Integer choicenbmult;///< user's choice for number of block multiplications (<0 -> automatic determination, min 6)
  Integer nblockmult;  ///< number of blockmultiplications in one restart
  Integer nlanczvecs;  ///< number of columns of storage matrix X carrying "meaningful" Ritz vectors 
  Integer retlanvecs;  ///< user's choice for number of returend Ritz vectors (<0 -> nlanczvecs)
  
  //--- global variables to allow reastarting
  Integer neigfound;  ///< number of eigenvalues known, these are in the first neigfound columns of X
  Integer blocksz;    ///< working block of lanczos vectors X(:,neigfound:neigfound+blocksz-1)
  Integer iter;       ///< iteration counter for Lanczos restarts (interval guess+computation)
  Integer nmult;      ///< number of single vector multiplications with matrix
  
  Real errc;   ///< error accumulation
  Real eps;    ///< relative precision
  Real mcheps; ///< machine precision (computed in constructor)
  Real maxval; ///< current approximation of largest eigenvalue
  Real minval; ///< current approximation of smallest eigenvalue
  Real polval; ///< the Chebychef polynomial will have this value at maxval
  
  Matrix X; ///< provides storage for the Lanczos vectors during computation
  Matrix C; ///< X^tAX, intermediate eigenvalue computations, orthogonalizations, etc
  Matrix d; ///< diagonals (Ritz values)
  Matrix e; ///< error bounds
  Matrix Xqr; ///< for complete orthogonalization with Householder QR
  Matrix u; ///< temporary matrix
  Matrix v; ///< temporary matrix 
  Matrix w; ///< temporary matrix
  Matrix minvec;  ///< temporary matrix for guessing minimal eigenvalue
  
  const Lanczosmatrix* bigmatrix; ///< pointer giving the (virtual) input matrix
    
  int stop_above;    ///< 1 if algorithm is to stop after upper bound is exceeded
  Real upper_bound;  ///< stop if current maximum Ritz value exceeds this value
  
  Integer ncalls;    ///< counts number of calls to Lanczpol
  Integer mymaxj;    ///< maximum amount of storage (columns) provided in X and C

  CH_Tools::GB_rand randgen;  ///< local random number generator 

  CH_Tools::Clock myclock;           ///< for time measurements 
  CH_Tools::Microseconds time_mult;  ///< for each restart, the time spent in lanczosmult
  CH_Tools::Microseconds time_mult_sum;  ///< sum over time_mult for all restarts
  CH_Tools::Microseconds time_iter; ///< time spent in last lanczos iteration
  CH_Tools::Microseconds time_sum; ///< time spent in this call to compute()

  int print_level;   ///<  level of iteration information that should be displayed
  std::ostream* myout;   ///< everything is output to *myout, default: myout=&cout (may be 0 for no output)

  /// output of error messages
  void error(const char *s){
    if (myout) 
      (*myout)<<"Lanczos Error: "<<s<<std::endl;
  }

  /// compute a guess for minimum and maximum eigenvalue
  int guess_extremes(Integer nproposed);
    
  /// main routine performing the Lancozs iterations
  int dhua(Integer nreig,Integer nproposed, Integer maxmult);
    
  /// block lanzos multiplication without Chebycheff spectral transformation
  int bklanc(Integer neigfound,Integer blocksz,Integer s,
               const Matrix &d,Matrix &C,Matrix &X,
               Matrix &e,Matrix &u,Matrix &v);

  /// block lanzos multiplication, with complete QR orthogonalization, without Chebycheff spectral transformation
  int bkqrlanc(Integer neigfound,Integer blocksz,Integer s,
               const Matrix &d,Matrix &C,Matrix &X,
               Matrix &e,Matrix &u,Matrix &v);

  /// block lanzos multiplication with Chebycheff spectral transformation
  int bklanccheb(Integer neigfound,Integer blocksz,Integer s,
                    Matrix &C,Matrix &X,Matrix &e,Matrix &v);

  /// block lanzos multiplication with complete QR orthogonalization and Chebycheff spectral transformation
  int bkqrlanccheb(Integer neigfound,Integer blocksz,Integer s,
                    Matrix &C,Matrix &X,Matrix &e,Matrix &v);
  
  int pch(Integer q,Integer neigfound,Integer nreig,Integer nconv,
	  Integer &blocksz,Integer &s,Integer iter,Integer sbs,const Matrix& d);

  /// compute norms of deviations of the Ritz vectors from being eigenvectors
  int err(Integer neigfound,Integer blocksz,const Matrix& X,Matrix &e);

  /// check convergence of maximum Ritz value/vector
  int cnvtst(Integer neigfound,Integer blocksz,Real& errc,Real eps,
               const Matrix& d,const Matrix &e,Integer &nconv);
    
  /// compute eigenvalues of current (block) tridiagonalization
  int eigen(Integer neigfound,Integer blocksz,Integer sbsz,Matrix& C,
	    Matrix& d,Matrix& u,Matrix& v,Real& af);

  /// compute matrix for eigenvalue computation into C
  int sectn(Matrix& X,Integer neigfound,Integer blocksz,
	    Matrix& C,Matrix& d,Matrix& u,Matrix& v,
	    Real& af);

  /// rotate the eigenvectors of the largest and smallest eigenvalues of the tridiagonal matrix.
  int rotate_extremes(Integer neigfound,Integer sbs,Matrix& d,
		      const Matrix& C,Matrix &X,Matrix& v);
  
  /// rotate the lanzos vectors to Ritz vectors
  int rotate(Integer neigfound,Integer sbs,Integer l,
	     const Matrix& C,Matrix &X,Matrix& v);

  /// assign a random vector to column j of X
  int random(Matrix& X,Integer j);

  /// orthonormalize columns X(:,offset:offset+l_blocksz-1) to all previous columns
  Integer orthog(Integer offset,Integer blocksz,Matrix &X,Matrix& B);

  /// apply spectral transformation by using a Chebycheff polynomial on the block multiplications
  int blockcheby(Integer col_offset,const Matrix& X,Matrix& v);

  /// compute the same polynomial as in blockcheby but for the scalar value xval
  Real scalarcheby(Real xval);

  /// tridiagonalize a symmetric blockdiagonal matrix
  int tred2(Integer n,const Matrix& C,Matrix& u,Matrix& v,Matrix& Z);

  /// compute the eigenvalues of a tridiagonal matrix 
  int tql2(Integer n,Matrix &u,Matrix &v, Matrix& Z);

public:
  /// intialize all to default values
  Lanczpol();
  /// destructor, nothing particular
  ~Lanczpol();


  /** @name Set and Get Parameters
       
      There should be no need to set any parameters, default values should be
      available and reasonable.
  */
  //@{
  
  /// set a guess on the value of the smallest   
  void set_mineig(Real ie){minval=ie;}
  /// set an upper bound on the number of matrix vector multiplications
  void set_maxmult(Integer mop){maxop=mop;}
  /// set an upper bpound on the number of restarts
  void set_maxiter(Integer mi){maxiter=mi;}
  /// set relative precision requirement for termination
  void set_relprec(Real relprec){eps=relprec;}  
  /// set the degree of the Chebycheff polynomial for the spectral transformation
  void set_nchebit(Integer cheb){choicencheb=cheb;} 
  /// set maximum number of block multiplications within one restart
  void set_nblockmult(Integer nb){choicenbmult=nb;}
  
  /// allow the algorithm to stop as soon as the maximum Ritz value exceeds the value ub
  void enable_stop_above(Real ub){stop_above=1; upper_bound=ub;}
  /// do not allow premature termination as in enable_stop_above()
  void disable_stop_above(){stop_above=0;}
  
  /// set an upper bound on the number of vectors returned in get_lanczosvecs()
  void set_retlanvecs(Integer nl){retlanvecs=nl;} 
  
  /// returns the Lanczos vectors of the last call with their Ritz values
  int get_lanczosvecs(Matrix& val,Matrix& vecs) const;
  
  /// returns current relative precision requirement
  Real get_relprec(void){return eps;} //relative precision
  
  /// returns the error code of the last call
  int get_err() const{return ierr;}
  /// returns the number of restarts of the last call
  Integer get_iter() const{return iter;}
    /// returns the number of matrix-vector multiplications of the last call
  Integer get_nmult() const{return nmult;}

    //@}
    
    /// the main routine: compute the nreig maximum eigenvalues of the matrix specified by bigmat
    int compute(const Lanczosmatrix* bigmat, ///< the symmetric matrix
		Matrix& eigval,        ///< on output: converged eigenvalues
		Matrix& eigvec,        ///< on output: eigenvectors to eigval, on input (optional): starting vectors
                Integer nreig,         ///< number of maximal eigenvalues to be computed
		Integer in_blocksz=0,  ///< size of a block, if block Lanczos is used
		Integer maxcol=0       ///< maximum number of columns that may be used
		);



  /** @name Input/Output
      
  */
  //@{
  
  /// set output stream and level of detail of log output (for debugging) 
  void set_out(std::ostream* o=0,int pril=1)
  {myout=o; print_level=pril;}


  /// save all data in out so that the current state can be recovered completely by restore()    
  std::ostream& save(std::ostream& out) const;
  
  /// restore the data from in where it was stored by save()
  std::istream& restore(std::istream& in);
  //@}
};

  //@}

}

#endif

