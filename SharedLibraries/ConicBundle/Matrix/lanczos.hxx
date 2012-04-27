/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/lanczos.hxx

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



#ifndef CH_MATRIX_CLASSES__LANCZOS_HXX
#define CH_MATRIX_CLASSES__LANCZOS_HXX

/**  @file lanczos.hxx
    @brief Header declaring the abstract classes CH_Matrix_Classes::Lanczosmatrix and CH_Matrix_Classes::Lanczos 
    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg

*/

#ifndef CH_MATRIX_CLASSES__MATRIX_HXX
#include "matrix.hxx"
#endif

namespace CH_Matrix_Classes {

/**@defgroup Lanczos_Interface Lanczos Interface and Classes 
   @brief Routines for computing a few extremal eigenvalues of large structured real symmetric matrices
*/
  //@{

  /** @brief Abstract base class for supplying the input matrix for Lanzcosmethods.
 
      It provides a virtual symmetric matrix with a routine for a 
      multiplying it with a matrix and a guess on the
      number of flops involved in this multiplication.
 
  */
  
  class Lanczosmatrix
  {
  public:
    virtual ~Lanczosmatrix(){}

    ///returns the order of the (virtual) symmetric matrix
    virtual Integer lanczosdim() const=0;
  
    ///returns a rough estimate on the number of flops needed by lanczosmult() for a  vector
    virtual Integer lanczosflops() const=0;

    ///computes  B = (*this) * A; A and B must not be the same object!
    virtual int lanczosmult(const Matrix& A,Matrix& B) const=0;
  };


  /** @brief Abstract interface to Lanzcos methods for computing a few extremal eigenvalues given via a Lanczosmatrix

  */

  class Lanczos
  {
  public:
    ///
    virtual ~Lanczos(){}
    
   /** @name Set and Get Parameters
       
       There should be no need to set any parameters, default values should be
       available and reasonable.
   */
    //@{
  
    /// set a guess on the value of the smallest 
    virtual void set_mineig(Real ie)=0;
    /// set an upper bound on the number of matrix vector multiplications
    virtual void set_maxmult(Integer mop)=0;
    /// set an upper bpound on the number of restarts
    virtual void set_maxiter(Integer mi)=0;
    /// set relative precision requirement for termination
    virtual void set_relprec(Real relprec)=0; 
    /// set maximum number of block multiplications within one restart
    virtual void set_nblockmult(Integer nb)=0;
    /// set the degree of the Chebycheff polynomial for the spectral transformation
    virtual void set_nchebit(Integer nc)=0;
    /// allow the algorithm to stop as soon as the maximum Ritz value exceeds the value ub
    virtual void enable_stop_above(Real ub)=0;
    /// do not allow premature termination as in enable_stop_above()
    virtual void disable_stop_above()=0;

    /// set an upper bound on the number of vectors returned in get_lanczosvecs()
    virtual void set_retlanvecs(Integer nl)=0; 

    /// returns the first neigfound+blocksz Lanczos vectors and their Ritz values of the last call
    virtual int get_lanczosvecs(Matrix& val,Matrix& vecs) const=0;
    
    /// returns current relative precision requirement
    virtual Real get_relprec(void)=0;

    /// returns the error code of the last call
    virtual int get_err() const=0;
    
    /// returns the number of restarts of the last call
    virtual Integer get_iter() const=0;
    
    /// returns the number of matrix-vector multiplications of the last call
    virtual Integer get_nmult() const=0;

    //@}
    
    /// compute the nreig maximum eigenvalues of the matrix specified by bigmat
    virtual int compute(const Lanczosmatrix* bigmat, ///< the symmetric matrix
                        Matrix& eigval,        ///< on output: converged eigenvalues
                        Matrix& eigvec,        ///< on output: eigenvectors to eigval, on input (optional): starting vectors
			Integer nreig,         ///< number of maximal eigenvalues to be computed
                        Integer in_blocksz=0,  ///< size of a block, if block Lanczos is used
                        Integer maxcol=0       ///< maximum number of columns that may be used
			)=0;

   /** @name Input/Output
       
   */
    //@{
  
    /// set output stream and level of detail of log output (for debugging) 
    virtual void set_out(std::ostream* out=0,int print_level=1)=0;

    /// save all data in out so that the current state can be recovered completely by restore()    
    virtual std::ostream& save(std::ostream& out) const =0;

    /// restore the data from in where it was stored by save()
    virtual std::istream& restore(std::istream& in) =0;

    //@}
  };

  //@}
  
}

#endif

