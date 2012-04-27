/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/indexmat.hxx

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



#ifndef CH_MATRIX_CLASSES__INDEXMAT_HXX
#define CH_MATRIX_CLASSES__INDEXMAT_HXX

/**  @file indexmat.hxx
    @brief Header declaring the classes CH_Matrix_Classes::Range and CH_Matrix_Classes::Indexmatrix for supporting integral matrices, typically needed for indexing purposes

    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg

*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#ifndef CH_MATRIX_CLASSES__MEMARRAY_HXX
#include "memarray.hxx"
#endif 
#ifndef CH_MATRIX_CLASSES__MYMATH_HXX
#include "mymath.hxx"
#endif


namespace CH_Matrix_Classes {


//everything involving a "Sparsesym" is implemented in "sparssym.cxx/hxx"
//everything else involving a "Sparsemat" is implemented in "sparsmat.cxx/hxx"
//everything else involving a "Symmatrix" is implemented in "symmat.cxx/hxx"
//everything else involving a "Matrix" is implemented in "matrix.cxx/hxx"


/**@defgroup rangegroup  Range (of integer numbers with step size)
*/
  //@{

// **************************************************************************
//                               Range
// **************************************************************************

//specifies a range of numbers: {i: i=from+k*step for some k\in N_0 and i<=to}

/// allows to specify a range of integral values via (from, to, step) meaning {j=from+i*step:j in[from,to],i in {0,1,2,...}}
class Range
{
public:
  ///
  Integer from;
  ///
  Integer to;
  ///  
  Integer step;

  ///
  Range(Integer _from,Integer _to,Integer _step=1):from(_from),to(_to),step(_step){}
  ///
  ~Range(){} 
};

  //@}



// **************************************************************************
//                              Indexmatrix
// **************************************************************************

class Matrix;
class Symmatrix;
class Sparsemat;
class Sparsesym;

/**@defgroup Indexmatrixgroup Indexmatrix (dense, integer, m by n)
*/
  //@{

  /** @brief %Matrix class for integral values of type #Integer

      Internally a matrix of size nr x nc is stored in a one dimensional array of ::Integer variables,
      the elements are arranged in columnwise order 
      (a11,a21,...,anr1,a12,a22,...).

      Any matrix element can be indexed by (i,j), which internally refers to m[i+j*nr], 
      or directly by the one dimensional index (i+j*nr). The latter view directly
      correspond to the vec() operator often used in the linear algebra literature,
      i.e., the matrix is transformed to a vector by stacking the columns on top
      of each other.
   */
class Indexmatrix: protected Memarrayuser
{
  friend class Matrix;
  friend class Symmatrix;
  friend class Sparsemat;
  friend class Sparsesym;

        
private:
  static const Mtype mtype;    ///< used for MatrixError templates (runtime type information was not yet existing)
    Integer mem_dim;   ///< amount of memory currently allocated
    Integer nr,        ///< number of rows
            nc;        ///< number of columns
    Integer *m;        ///< pointer to store, order is columnwise (a11,a21,...,anr1,a12,a22,.....)

    bool is_init;      ///< flag whether memory is initialized, it is only used if CONICBUNDLE_DEBUG is defined

    /// initialize the matrix to a 0x0 matrix without storage
    inline void init_to_zero();   
    
public:

  //----------------------------------------------------
  //----  constructors, destructor, and initialization
  //----------------------------------------------------

  /** @name Constructors, Destructor, and Initialization (Members)
   */
  //@{

  /// empty matrix
  inline Indexmatrix();                  

  /// copy constructor, *this=d*A
  inline Indexmatrix(const Indexmatrix& A,Integer d=1);

  /// generate a column vector holding the indices of this #Range
  inline Indexmatrix(const Range&);           

  /** @brief generate a matrix of size nr x nc but WITHOUT initializing the memory
      
      If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use 
      set_init() via matrix.set_init(true) in order to avoid warnings concerning improper 
      initialization
  */ 
  inline Indexmatrix(Integer nr,Integer nc);         

  /// generate a matrix of size nr x nc initializing all elements to the value d
  inline Indexmatrix(Integer nr,Integer nc,Integer d);  

  /// generate a matrix of size nr x nc initializing the elements from the (one dimensional) array dp with increment incr
  inline Indexmatrix(Integer nr,Integer nc,const Integer* dp,Integer incr=1); 

  /// copy a std::vector<#Integer> to a column vector of the same size and content
  inline Indexmatrix(const std::vector<Integer>& vec);

  ///
  ~Indexmatrix(){memarray->free(m);}
  
#if (CONICBUNDLE_DEBUG>=1)
  /// after external initialization, call matrix.set_init(true) (not needed if CONICBUNDLE_DEBUG is undefined)
  void set_init(bool i){is_init=i;}
  /// returns true if the matrix has been declared initialized (not needed if CONICBUNDLE_DEBUG is undefined)
  bool get_init() const {return is_init;}
#else 
  /// after external initialization, call matrix.set_init(true) (not needed if CONICBUNDLE_DEBUG is undefined) 
  void set_init(bool /* i */){}
  /// returns true if the matrix has been declared initialized (not needed if CONICBUNDLE_DEBUG is undefined)
  bool get_init() const {return true;}
#endif
  
  /// initialize to *this=A*d
  inline Indexmatrix& init(const Indexmatrix &A,Integer d=1);    

  /// initialize *this to a column vector holding the indices of #Range
  Indexmatrix& init(const Range&);

  /// intialize *this to a matrix of size nr x nc initializing all elements to the value d
  inline Indexmatrix& init(Integer nr,Integer nc,Integer d);                    

  /// generate a matrix of size nr x nc initializing the elements from the (one dimensional) array dp with increment incr
  inline Indexmatrix& init(Integer nr,Integer nc,const Integer *dp,Integer incr=1);

  /// use std::vector<#Integer> to initialize this to a column vector of the same size and content
  inline Indexmatrix& init(const std::vector<Integer>& vec);
    
  /** @brief resize the matrix to nr x nc elements but WITHOUT initializing the memory
      
      If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use 
      set_init() via matrix.set_init(true) in order to avoid warnings concerning improper 
      initialization
  */ 
  void newsize(Integer nr,Integer nc);  //resize matrix without initialization

  //@}

  /** @name Conversions from other Matrix Classes (Members)
   */
  //@{

  /// copy with rounding
  Indexmatrix(const Matrix&);          
  ///copy with rounding
  Indexmatrix(const Symmatrix&);      
  ///copy with rounding
  Indexmatrix(const Sparsemat&);      
  ///copy with rounding
  Indexmatrix(const Sparsesym&);      

  //@}

  //----------------------------------------------------
  //----  size and type information
  //----------------------------------------------------

  /** @name Size and Type Information (Members)
   */
  //@{

  /// returns the number of rows in _nr and the number of columns in _nc
  void dim(Integer& _nr, Integer& _nc) const {_nr=nr; _nc=nc;}

  /// returns the dimension rows * columns when the matrix is regarded as a vector
  Integer dim() const {return nr*nc;}

  /// returns the row dimension
  Integer rowdim() const {return nr;}

  /// returns the column dimension
  Integer coldim() const {return nc;}

  /// returns the type of the matrix, MTindexmatrix
  Mtype get_mtype() const {return mtype;}

  //@}

  //--------------------------------
  //----  Indexing and Submatrices
  //--------------------------------


  /** @name Indexing and Submatrices (Members)
   */
  //@{
  
  /// returns reference to element (i,j) of the matrix (rowindex i, columnindex j)
  inline Integer& operator()(Integer i,Integer j);      

  /// returns reference to element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
  inline Integer& operator()(Integer i);      //index to vector of stacked columns

  /// returns value of element (i,j) of the matrix (rowindex i, columnindex j)
  inline Integer operator()(Integer i,Integer j) const; 

  /// returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
  inline Integer operator()(Integer i) const; //index to vector of stacked columns

  /// returns a new submatrix as indexed by vecrow and veccol, A(i,j)=(*this)(vecrow(i),veccol(j)) for 0<=i<vecrow.dim(), 0<=j<veccol.dim() 
  Indexmatrix operator()(const Indexmatrix& vecrow,const Indexmatrix& veccol) const;

  /// returns a new matrix B of the same shape as A with B(i,j)=(*this)(A(i),A(j)) for 0<=i<A.rowdim(), 0<=j<A.coldim() 
    Indexmatrix operator()(const Indexmatrix& A) const;
      //*this matrix is interpreted as vector, the resulting matrix has
      //has the same shape as A
  
  /// returns reference to element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
  inline Integer& operator[](Integer i);      //{return (*this)(i);}

  /// returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
  inline Integer operator[](Integer i) const; //{return (*this)(i);}
  
  /// returns column i copied to a new matrix 
  Indexmatrix col(Integer i) const;          //returns column i as column vector
  /// returns row i copied to a new matrix 
  Indexmatrix row(Integer i) const;          //returns row i as row vector
  /// returns a matrix of size this->rowdim() x vec.dim(), with column i a copy of column vec(i) of *this
  Indexmatrix cols(const Indexmatrix &vec) const;  //returns cols as indexed by vec
  
  /// returns a matrix of size vec.dim() x this->rowdim(), with row i a copy of row vec(i) of *this
  Indexmatrix rows(const Indexmatrix &vec) const;  //returns rows as indexed by vec


  /// keeps upper triangle starting with diagonal d; set (*this)(i,j)=0 for 0<=i<row dimension, 0<=j<min(i+d,column dimension), returns *this
  Indexmatrix& triu(Integer d=0);    // (*this)(i,j)=0 for i<j+d

  /// keeps lower triangle starting with diagonal d; set (*this)(i,j)=0 for 0<=i<row dimension, max(0,i+d)<=j<column dimension, returns *this
  Indexmatrix& tril(Integer d=0);    // (*this)(i,j)=0 for i>j+d

  /// assigns A to a submatrix of *this,  (*this)(vecrow(i),veccol(j))=A(i,j) for 0<=i<vecrow.dim(), 0<=j<veccol.dim()
  Indexmatrix& subassign(const Indexmatrix& vecrow,const Indexmatrix& veccol,
		    const Indexmatrix& A);
     //(*this)(vecrow(i),veccol(j))=A(i,j) for all i,j

  /// assigns vector A to a subvector of *this,  (*this)(vec(i))=A(i) for 0<=i<vec.dim(), *this, vec, and A may be rectangular matrices, their dimesions are not changed, returns *this
  Indexmatrix& subassign(const Indexmatrix& vec,const Indexmatrix& A);

  /// all rows indexed by vector ind are deleted, no row should appear twice in ind, remaining rows are moved up keeping their order, returns *this
  Indexmatrix& delete_rows(const Indexmatrix& ind);
  /// all colmuns indexed by vector ind are deleted, no column should appear twice in ind, remaining columns are moved up keeping their order, returns *this
  Indexmatrix& delete_cols(const Indexmatrix& ind);

  /// insert the row vector v before row i, 0<=i<= row dimension, for i==row dimension the row is appended below; appending to a 0x0 matrix is allowed, returns *this
  Indexmatrix& insert_row(Integer i,const Indexmatrix& v); // 0<=i<=nr, v vector
  /// insert a column before column i, 0<=i<= column dimension, for i==column dimension the column is appended at the right; appending to a 0x0 matrix is allowed, returns *this
  Indexmatrix& insert_col(Integer i,const Indexmatrix& v); // 0<=i<=nc, v vector

  /// (*this) is set to a column vector of length min{max{0,n},dim()}; usually used to truncate a vector, returns *this
  inline Indexmatrix& reduce_length(Integer n);  

  /// concats matrix A to the right of *this, A or *this may be the 0x0 matrix initally, returns *this
  Indexmatrix& concat_right(const Indexmatrix& A);
  /// concats matrix A to the bottom of *this, A or *this may be the 0x0 matrix initally, returns *this
  Indexmatrix& concat_below(const Indexmatrix& A);
  /// concat value d at the bottom of *this, *this must be a column vector or the 0x0 matrix, returns *this 
  Indexmatrix& concat_below(Integer d);   //only for column vectors (nr==1)
  /// concat value d at the right of *this, *this must be a row vector or the 0x0 matrix, returns *this 
  Indexmatrix& concat_right(Integer d);   //only for row vectors (nc==1)


  /// returns the current address of the internal value array; use cautiously, do not use delete!
  Integer* get_store() {return m;}       //use cautiously, do not use delete!
  /// returns the current address of the internal value array; use cautiously!
  const Integer* get_store() const {return m;}   //use cautiously

  //@}


  /** @name Indexing and Submatrices (Friends)
   */
  //@{


  /// returns a column vector v consisting of the elements v(i)=(*this)(i,i), 0<=i<min(row dimension,column dimension) 
  friend Indexmatrix diag(const Indexmatrix& A);      //=(A(1,1),A(2,2),...)^t

  /// retuns a matrix that keeps the upper triangle of A starting with diagonal d, i.e., (i,j)=A(i,j) for 0<=i<row dimension, max(0,i+d)<=j<column dimension, and sets (i,j)=0 otherwise
  friend inline Indexmatrix triu(const Indexmatrix& A,Integer d=0);

  /// retuns a matrix that keeps the lower triangle of A starting with diagonal d, i.e., (i,j)=A(i,j) for 0<=i<row dimension, 0<=j<min(i+d+1,column dimension), and sets (i,j)=0 otherwise
  friend inline Indexmatrix tril(const Indexmatrix& A,Integer d=0);

  /// returns a new matrix [A, B], i.e., it concats matrices A and B rowwise; A or B may be a 0x0 matrix
  friend inline Indexmatrix concat_right(const Indexmatrix& A,const Indexmatrix& B);
  /// returns the matrix [A; B], i.e., it concats matrices A and B columnwise; A or B may be a 0x0 matrix
  friend inline Indexmatrix concat_below(const Indexmatrix& A,const Indexmatrix& B);

  /// swap the content of the two matrices A and B (involves no copying)
  friend void swap(Indexmatrix &A,Indexmatrix &B);

  //@}

  
  //------------------------------
  //----  BLAS-like Routines
  //------------------------------

  /** @name BLAS-like Routines (Members)
   */
  //@{

  ///sets *this=d*A and returns *this
  Indexmatrix& xeya(const Indexmatrix& A,Integer d=1); 
  ///sets *this+=d*A and returns *this  
  Indexmatrix& xpeya(const Indexmatrix& A,Integer d=1);
    
  //@}

  /** @name BLAS-like Routines (Friends)
   */
  //@{

  ///returns x= alpha*y+beta*x, where y may be transposed (ytrans=1); if beta==0. then x is initialized to the correct size
  friend Indexmatrix& xbpeya(Indexmatrix& x,const Indexmatrix& y,Integer alpha=1,Integer beta=0,int ytrans=0);
  
  ///returns x= alpha*y+beta*z; x is initialized to the correct size
  friend Indexmatrix& xeyapzb(Indexmatrix& x,const Indexmatrix& y,const Indexmatrix& z,Integer alpha=1,Integer beta=1);
  
  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0 then C is initialized to the correct size
  friend Indexmatrix& genmult(const Indexmatrix& A,const Indexmatrix& B,Indexmatrix& C,
			 Integer alpha=1,Integer beta=0,int atrans=0,int btrans=0);
    
  //@}
    

  //-----------------------------------
  //----  usual arithmetic operators
  //-----------------------------------

  /** @name Usual Arithmetic Operators (Members)
   */
  //@{

  ///
  inline Indexmatrix& operator=(const Indexmatrix &A);
  ///
  inline Indexmatrix& operator*=(const Indexmatrix &s);
  ///
  inline Indexmatrix& operator+=(const Indexmatrix &v);
  ///
  inline Indexmatrix& operator-=(const Indexmatrix &v);
  /// ATTENTION: this is redefined as the Hadamard product, (*this)(i,j)=(*this)(i,j)*A(i,j) for all i,j
  inline Indexmatrix& operator%=(const Indexmatrix &A);   //Hadamard product!!
  ///
  inline Indexmatrix operator-() const;
  
  /// 
  inline Indexmatrix& operator*=(Integer d);
  /// ATTENTION: d is NOT checked for 0
  inline Indexmatrix& operator/=(Integer d);
  /// sets (*this)(i,j)\%=d for all i,j in the modulo meaning of \%, ATTENTION: d is NOT checked for 0
  inline Indexmatrix& operator%=(Integer d); 
  /// sets (*this)(i,j)+=d for all i,j
  inline Indexmatrix& operator+=(Integer d);
  /// sets (*this)(i,j)-=d for all i,j
  inline Indexmatrix& operator-=(Integer d);
  
  ///transposes itself (cheap for vectors, expensive for matrices)
  Indexmatrix& transpose();
  
  //@}

  /** @name Usual Arithmetic Operators (Friends)
   */
  //@{

  ///
  friend inline Indexmatrix operator*(const Indexmatrix &A,const Indexmatrix& B);
  ///
  friend inline Indexmatrix operator+(const Indexmatrix &A,const Indexmatrix& B);
  ///
  friend inline Indexmatrix operator-(const Indexmatrix &A,const Indexmatrix& B);
  /// ATTENTION: this is redefined as the Hadamard product, (*this)(i,j)=(*this)(i,j)*A(i,j) for all i,j
  friend inline Indexmatrix operator%(const Indexmatrix &A,const Indexmatrix& B);
  /// 
  friend inline Indexmatrix operator*(const Indexmatrix &A,Integer d);
  ///
  friend inline Indexmatrix operator*(Integer d,const Indexmatrix &A);
  /// ATTENTION: d is NOT checked for 0
  friend inline Indexmatrix operator/(const Indexmatrix& A,Integer d);
  /// sets (i,j)=A(i,j)\%d for all i,j in the modulo meaning of \%, ATTENTION: d is NOT checked for 0
  friend inline Indexmatrix operator%(const Indexmatrix& A,Integer d);
  
  /// returns (i,j)=A(i,j)+d for all i,j
  friend inline Indexmatrix operator+(const Indexmatrix& A,Integer d);
  /// returns (i,j)=A(i,j)+d for all i,j
  friend inline Indexmatrix operator+(Integer d,const Indexmatrix& A);
  /// returns (i,j)=A(i,j)-d for all i,j
  friend inline Indexmatrix operator-(const Indexmatrix& A,Integer d);
  /// returns (i,j)=d-A(i,j) for all i,j
  friend inline Indexmatrix operator-(Integer d,const Indexmatrix& A);

  ///
  friend Indexmatrix transpose(const Indexmatrix& A);

  //@}

  //------------------------------------------
  //----  Connections to other Matrix Classes
  //------------------------------------------

  /** @name Connections to other Classes (Friends)
   */
  //@{

  /// interpret A as a vector and copy it to a std::vector<int> which is also returned 
  friend inline std::vector<int>& assign(std::vector<int>& vec,
					    const Indexmatrix& A); 
  /// interpret A as a vector and copy it to a std::vector<long> which is also returned 
  friend inline std::vector<long>& assign(std::vector<long>& vec,
					     const Indexmatrix& A); 

  //@}


  //------------------------------
  //----  Elementwise Operations
  //------------------------------

  /** @name Elementwise Operations (Members)
   */
  //@{

  /// resize *this to an nr x nc matrix and assign to (i,j) a random number uniformly from [lowerb,upperb] for all i,j
  Indexmatrix& rand(Integer nr,Integer nc,Integer lowerb,Integer upperb,CH_Tools::GB_rand* random_generator=0);
  /// shuffle the elements randomly (does not change dimensions)
  Indexmatrix& shuffle(CH_Tools::GB_rand* random_generator=0);
  /// using ::sign assign (*this)(i,j)=sign((*this)(i,j)) for all i,j 
  Indexmatrix& sign(void);
  /// using ::abs assign (*this)(i,j)=abs((*this)(i,j)) for all i,j 
  Indexmatrix& abs(void);        

  //@}

  /** @name Elementwise Operations (Friends)
   */
  //@{

  /// return a nr x nc matrix with (i,j) assigned a random number uniformly from [lowerb,upperb] for all i,j
  friend inline Indexmatrix rand(Integer nr,Integer nc,Integer lb,Integer ub,CH_Tools::GB_rand* random_generator=0);
  /// return a matrix of the same size as A with (i,j)=sign(A(i,j)) for all i,j, see also CH_Matrix_Classes::sign()
  friend inline Indexmatrix sign(const Indexmatrix& A);
  friend Indexmatrix abs(const Indexmatrix& A);                 

  //@}

  //----------------------------
  //----  Numerical Methods
  //----------------------------

  /** @name Numerical Methods (Members)
   */
  //@{

  ///scales each row i of (*this) by vec(i), i.e., (*this)=diag(vec)*(*this), and returns (*this)
  Indexmatrix& scale_rows(const Indexmatrix& vec);    
  ///scales each column i of (*this) by vec(i), i.e., (*this)=(*this)*diag(vec), and returns (*this)
  Indexmatrix& scale_cols(const Indexmatrix& vec); 

  //@}

  /** @name Numerical Methods (Friends)
   */
  //@{

  /// returns the sum of the diagonal elements A(i,i) over all i
  friend Integer trace(const Indexmatrix& A);               //=sum(diag(A))
  ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  friend inline Integer ip(const Indexmatrix& A, const Indexmatrix& B); //=trace(B^t*A)
  ///returns the Frobenius norm of A, i.e., the square root of the sum of A(i,j)*A(i,j) over all i,j 
  friend inline Real norm2(const Indexmatrix& A); 

  
  ///returns a row vector holding the sum over all rows, i.e., (1 1 ... 1)*A
  friend Indexmatrix sumrows(const Indexmatrix& A);   
  ///returns a column vector holding the sum over all columns, i.e., A*(1 1 ... 1)^T
  friend Indexmatrix sumcols(const Indexmatrix& A);
  ///returns the sum over all elements of A, i.e., (1 1 ... 1)*A*(1 1 ... 1)^T
  friend Integer sum(const Indexmatrix& A);     

  //@}


  //---------------------------------------------
  //----  Comparisons / Max / Min / sort / find
  //---------------------------------------------

  /** @name Find (Members)
   */
  //@{

  /// returns an Indexmatrix ind so that (*this)(ind(i)) 0<=i<ind.dim() runs through all nonzero elements  
  Indexmatrix find() const;   //finds nonzeros
  /// returns an Indexmatrix ind so that (*this)(ind(i)) 0<=i<ind.dim() runs through all elements having value num 
  Indexmatrix find_number(Integer num=0) const; 

  //@}

  /** @name Comparisons, Max, Min, Sort, Find (Friends)
   */
  //@{

  /// returns a matrix having elements (i,j)=Integer(A(i,j)<B(i,j)) for all i,j
  friend Indexmatrix operator<(const Indexmatrix &A,const Indexmatrix &B);
  /// returns a matrix having elements (i,j)=Integer(A(i,j)>B(i,j)) for all i,j
  friend inline Indexmatrix operator>(const Indexmatrix &A,const Indexmatrix &B);
  /// returns a matrix having elements (i,j)=Integer(A(i,j)<=B(i,j)) for all i,j
  friend Indexmatrix operator<=(const Indexmatrix &A,const Indexmatrix &B);
  /// returns a matrix having elements (i,j)=Integer(A(i,j)>=B(i,j)) for all i,j
  friend inline Indexmatrix operator>=(const Indexmatrix &A,const Indexmatrix &B);
  /// returns a matrix having elements (i,j)=Integer(A(i,j)==B(i,j)) for all i,j
  friend Indexmatrix operator==(const Indexmatrix &A,const Indexmatrix &B);
  /// returns a matrix having elements (i,j)=Integer(A(i,j)!=B(i,j)) for all i,j
  friend Indexmatrix operator!=(const Indexmatrix &A,const Indexmatrix &B);
  /// returns a matrix having elements (i,j)=Integer(A(i,j)<d) for all i,j
  friend Indexmatrix operator<(const Indexmatrix &A,Integer d);
  /// returns a matrix having elements (i,j)=Integer(A(i,j)>d) for all i,j
  friend Indexmatrix operator>(const Indexmatrix &A,Integer d);
  /// returns a matrix having elements (i,j)=Integer(A(i,j)<=d) for all i,j
  friend Indexmatrix operator<=(const Indexmatrix &A,Integer d);
  /// returns a matrix having elements (i,j)=Integer(A(i,j)>=d) for all i,j
  friend Indexmatrix operator>=(const Indexmatrix &A,Integer d);
  /// returns a matrix having elements (i,j)=Integer(A(i,j)==d) for all i,j
  friend Indexmatrix operator==(const Indexmatrix &A,Integer d);
  /// returns a matrix having elements (i,j)=Integer(A(i,j)!=d) for all i,j
  friend Indexmatrix operator!=(const Indexmatrix &A,Integer d);
  /// returns a matrix having elements (i,j)=Integer(d<A(i,j)) for all i,j
  friend inline Indexmatrix operator<(Integer d,const Indexmatrix &A);
  /// returns a matrix having elements (i,j)=Integer(d>A(i,j)) for all i,j
  friend inline Indexmatrix operator>(Integer d,const Indexmatrix &A);
  /// returns a matrix having elements (i,j)=Integer(d<=A(i,j)) for all i,j
  friend inline Indexmatrix operator<=(Integer d,const Indexmatrix &A);
  /// returns a matrix having elements (i,j)=Integer(d>=A(i,j)) for all i,j
  friend inline Indexmatrix operator>=(Integer d,const Indexmatrix &A);
  /// returns a matrix having elements (i,j)=Integer(d==A(i,j)) for all i,j
  friend inline Indexmatrix operator==(Integer d,const Indexmatrix &A);
  /// returns a matrix having elements (i,j)=Integer(d!=A(i,j)) for all i,j
  friend inline Indexmatrix operator!=(Integer d,const Indexmatrix &A);
   
  /// returns a row vector holding in each column the minimum over all rows in this column
  friend Indexmatrix minrows(const Indexmatrix& A);   
  /// returns a column vector holding in each row the minimum over all columns in this row
  friend Indexmatrix mincols(const Indexmatrix& A);   
  /// returns the minimum value over all elements of the matrix
  friend Integer min(const Indexmatrix& A,Integer *iindex=0,Integer *jindex=0);
  /// returns a row vector holding in each column the maximum over all rows in this column
  friend Indexmatrix maxrows(const Indexmatrix& A);   
  /// returns a column vector holding in each row the maximum over all columns in this row
  friend Indexmatrix maxcols(const Indexmatrix& A);   
  /// returns the maximum value over all elements of the matrix
  friend Integer max(const Indexmatrix& A,Integer *iindex=0,Integer *jindex=0);
  
  /// returns an Indexmatrix ind so that vec(ind(0))<=vec(ind(1))<=...<=vec(ind(vec.dim()-1)) (vec may be rectangular)
  friend inline Indexmatrix sortindex(const Indexmatrix& vec);
  /// sets ind so that vec(ind(0))<=vec(ind(1))<=...<=vec(ind(vec.dim()-1)) (vec may be rectangular)
  friend void sortindex(const Indexmatrix& vec,Indexmatrix &ind);
  
  /// returns an Indexmatrix ind so that A(ind(i)) 0<=i<ind.dim() runs through all nonzero elements of A 
  friend inline Indexmatrix find(const Indexmatrix& A);
  /// returns an Indexmatrix ind so that A(ind(i)) 0<=i<ind.dim() runs through all elements of A having value num 
  friend inline Indexmatrix find_number(const Indexmatrix& A,Integer num=0);

  //@}

  //--------------------------------
  //----  Input / Output
  //--------------------------------

  /** @name Input, Output (Members)
   */
  //@{

  /** @brief displays a matrix in a pretty way for bounded screen widths; for variables of value zero default values are used.
   */
    void display(std::ostream& out, ///< output stream
		 int precision=0,   ///< not needed here, used for consistency with real valued matrices
		 int width=0,       ///< field width, default = precision+6
		 int screenwidth=0  ///< maximum number of characters in one output line, default = 80
		 ) const;

  //@}

  /** @name Input, Output (Friends)
   */
  //@{

  ///output format (all ::Integer values): nr nc \\n A(1,1) A(1,2) ... A(1,nc) \\n A(2,1) ... A(nr,nc) \\n
  friend std::ostream& operator<<(std::ostream& o,const Indexmatrix &A);
  ///input format (all ::Integer values): nr nc \\n A(1,1) A(1,2) ... A(1,nc) \\n A(2,1) ... A(nr,nc) \\n
  friend std::istream& operator>>(std::istream& i,Indexmatrix &A);
  
  //@}  

};

//@}

// **************************************************************************
//                make non inline friends available outside
// **************************************************************************

Indexmatrix diag(const Indexmatrix& A);      //=(A(1,1),A(2,2),...)^t
void swap(Indexmatrix &A,Indexmatrix &B);
Indexmatrix& xbpeya(Indexmatrix& x,const Indexmatrix& y,Integer alpha,Integer beta,int ytrans);
Indexmatrix& xeyapzb(Indexmatrix& x,const Indexmatrix& y,const Indexmatrix& z,Integer alpha,Integer beta);
Indexmatrix& genmult(const Indexmatrix& A,const Indexmatrix& B,Indexmatrix& C,
			 Integer alpha,Integer beta,int atrans,int btrans);
Indexmatrix transpose(const Indexmatrix& A);
Indexmatrix abs(const Indexmatrix& A);                 
Integer trace(const Indexmatrix& A);               //=sum(diag(A))
Indexmatrix sumrows(const Indexmatrix& A);   //=(1 1 1 ... 1)*A
Indexmatrix sumcols(const Indexmatrix& A);   //=A*(1 1 ... 1)^t
Integer sum(const Indexmatrix& A);         //=(1 1 ... 1)*A*(1 1 ... 1)^t
Indexmatrix operator<(const Indexmatrix &A,const Indexmatrix &B);
Indexmatrix operator<=(const Indexmatrix &A,const Indexmatrix &B);
Indexmatrix operator==(const Indexmatrix &A,const Indexmatrix &B);
Indexmatrix operator!=(const Indexmatrix &A,const Indexmatrix &B);
Indexmatrix operator<(const Indexmatrix &A,Integer d);
Indexmatrix operator>(const Indexmatrix &A,Integer d);
Indexmatrix operator<=(const Indexmatrix &A,Integer d);
Indexmatrix operator>=(const Indexmatrix &A,Integer d);
Indexmatrix operator==(const Indexmatrix &A,Integer d);
Indexmatrix operator!=(const Indexmatrix &A,Integer d);
Indexmatrix minrows(const Indexmatrix& A);   //min of each column (over the rows)
Indexmatrix mincols(const Indexmatrix& A);   //min of each row (over the columns)
Integer min(const Indexmatrix& A,Integer *iindex,Integer *jindex);
Indexmatrix maxrows(const Indexmatrix& A);   //similar
Indexmatrix maxcols(const Indexmatrix& A);   //similar
Integer max(const Indexmatrix& A,Integer *iindex,Integer *jindex);
Indexmatrix sortindex(const Indexmatrix& vec);
void sortindex(const Indexmatrix& vec,Indexmatrix &ind);
std::ostream& operator<<(std::ostream& o,const Indexmatrix &v);
std::istream& operator>>(std::istream& i,Indexmatrix &v);




// **************************************************************************
//                   implementation of inline functions
// **************************************************************************

inline void Indexmatrix::init_to_zero()
{
 nr=nc=0;mem_dim=0;m=0;
 chk_set_init(*this,1);
}

inline Indexmatrix& Indexmatrix::init(Integer inr,Integer inc,Integer d)
{
 newsize(inr,inc);
 mat_xea(nr*nc,m,d);
 chk_set_init(*this,1);
 return *this;
}

inline Indexmatrix& Indexmatrix::init(Integer inr,Integer inc,const Integer* d,Integer incr)
{
 newsize(inr,inc);
 if (incr==1) mat_xey(nr*nc,m,d);
 else mat_xey(nr*nc,m,Integer(1),d,incr);
 chk_set_init(*this,1);
 return *this;
}

inline Indexmatrix& Indexmatrix::init(const Indexmatrix& A,Integer d)
{
 return xeya(A,d);
}

inline Indexmatrix& Indexmatrix::init(const std::vector<Integer>& vec)
{
  newsize(Integer(vec.size()),1); chk_set_init(*this,1);
  for(Integer i=0;i<nr;i++) m[i]=vec[i];
  return *this;
}

inline Indexmatrix::Indexmatrix()
{
 init_to_zero();
}

inline Indexmatrix::Indexmatrix(const Indexmatrix &A,Integer d):Memarrayuser()
{
 init_to_zero();
 xeya(A,d);
}

inline Indexmatrix::Indexmatrix(Integer inr,Integer inc)
{
 init_to_zero();
 newsize(inr,inc);
}

inline Indexmatrix::Indexmatrix(const Range& range)
{
 init_to_zero();
 init(range);
}

inline Indexmatrix::Indexmatrix(Integer inr,Integer inc,Integer d)
{
 init_to_zero();
 init(inr,inc,d);
}

inline Indexmatrix::Indexmatrix(Integer inr,Integer inc,const Integer *d,Integer incr)
{
 init_to_zero();
 init(inr,inc,d,incr);
}

inline Indexmatrix::Indexmatrix(const std::vector<Integer>& vec)
{
 init_to_zero();
 init(vec);
}

inline Integer& Indexmatrix::operator()(Integer i,Integer j)
{
 chk_range(i,j,nr,nc);
 return m[j*nr+i];
}

inline Integer& Indexmatrix::operator()(Integer i) 
{
 chk_range(i,0,nr*nc,1);
 return m[i];
}

inline Integer Indexmatrix::operator()(Integer i,Integer j) const
{
 chk_range(i,j,nr,nc);
 return m[j*nr+i];
}

inline Integer Indexmatrix::operator()(Integer i) const
{
 chk_range(i,0,nr*nc,1);
 return m[i];
}

inline Integer& Indexmatrix::operator[](Integer i)
{return (*this)(i);}

inline Integer Indexmatrix::operator[](Integer i) const 
{return (*this)(i);}

inline Indexmatrix& Indexmatrix::reduce_length(Integer n) 
{ nr=min(nr*nc,max(Integer(0),n)); nc=1; return *this;}


inline Integer ip(const Indexmatrix& A, const Indexmatrix& B)
{
 chk_add(A,B);
 return mat_ip(A.nc*A.nr,A.m,B.m);
}

inline Real norm2(const Indexmatrix& A)
{
 chk_init(A);
 return ::sqrt(double(mat_ip(A.nc*A.nr,A.m,A.m)));
}

inline Indexmatrix& Indexmatrix::operator=(const Indexmatrix &A)
{ return xeya(A);}

inline Indexmatrix& Indexmatrix::operator*=(const Indexmatrix &A)
{ Indexmatrix C; return xeya(genmult(*this,A,C));}

inline Indexmatrix& Indexmatrix::operator+=(const Indexmatrix &A)
{ return xpeya(A); }

inline Indexmatrix& Indexmatrix::operator-=(const Indexmatrix &A)
{ return xpeya(A,-1); }

inline Indexmatrix& Indexmatrix::operator%=(const Indexmatrix &A)
{ chk_add(*this,A); mat_xhadey(nr*nc,m,A.m); return *this; }

inline Indexmatrix Indexmatrix::operator-() const
{ return Indexmatrix(*this,-1); }

inline Indexmatrix& Indexmatrix::operator*=(register Integer d)
{ chk_init(*this); mat_xmultea(nr*nc,m,d); return *this; }

inline Indexmatrix& Indexmatrix::operator/=(register Integer d)
{ chk_init(*this); mat_xdivea(nr*nc,m,d); return *this; }

inline Indexmatrix& Indexmatrix::operator%=(register Integer d)
{ chk_init(*this); mat_xmodea(nr*nc,m,d); return *this; }

inline Indexmatrix& Indexmatrix::operator+=(register Integer d)
{ chk_init(*this); mat_xpea(nr*nc,m,d); return *this; }

inline Indexmatrix& Indexmatrix::operator-=(register Integer d)
{ chk_init(*this); mat_xpea(nr*nc,m,-d); return *this; }


inline Indexmatrix operator*(const Indexmatrix &A,const Indexmatrix &B) 
    {Indexmatrix C; return genmult(A,B,C);}
inline Indexmatrix operator+(const Indexmatrix &A,const Indexmatrix &B)
    {Indexmatrix C; return xeyapzb(C,A,B,1,1);}
inline Indexmatrix operator-(const Indexmatrix &A,const Indexmatrix &B)
    {Indexmatrix C; return xeyapzb(C,A,B,1,-1);}
inline Indexmatrix operator%(const Indexmatrix &A,const Indexmatrix &B) 
    {Indexmatrix C(A); return C%=B;}
inline Indexmatrix operator*(const Indexmatrix &A,Integer d)          
    {return Indexmatrix(A,d);}
inline Indexmatrix operator*(Integer d,const Indexmatrix &A)
    {return Indexmatrix(A,d);}
inline Indexmatrix operator/(const Indexmatrix &A,Integer d)
    {Indexmatrix B(A); return B/=d;}
inline Indexmatrix operator%(const Indexmatrix &A,Integer d)
    {Indexmatrix B(A); return B%=d;}
inline Indexmatrix operator+(const Indexmatrix &A,Integer d)           
    {Indexmatrix B(A); return B+=d;}
inline Indexmatrix operator+(Integer d,const Indexmatrix &A)           
    {Indexmatrix B(A); return B+=d;}
inline Indexmatrix operator-(const Indexmatrix &A,Integer d)       
    {Indexmatrix B(A); return B-=d;}
inline Indexmatrix operator-(Integer d,const Indexmatrix &A)       
    {Indexmatrix B(A,-1); return B+=d;}

inline Indexmatrix triu(const Indexmatrix& A,Integer i)
          {Indexmatrix B(A); B.triu(i); return B;}
inline Indexmatrix tril(const Indexmatrix& A,Integer i)
          {Indexmatrix B(A); B.tril(i); return B;}

inline Indexmatrix concat_right(const Indexmatrix& A,const Indexmatrix& B)
          {Indexmatrix C(A.dim()+B.dim(),1);C=A;C.concat_right(B);return C;}
inline Indexmatrix concat_below(const Indexmatrix& A,const Indexmatrix& B)
          {Indexmatrix C(A.dim()+B.dim(),1);C=A;C.concat_below(B);return C;}
    
inline void swap(Indexmatrix &A,Indexmatrix &B)
{ 
 Integer *hm=A.m;A.m=B.m;B.m=hm;
 Integer hi=A.nr;A.nr=B.nr;B.nr=hi;
 hi=A.nc;A.nc=B.nc;B.nc=hi;
 hi=A.mem_dim;A.mem_dim=B.mem_dim;B.mem_dim=hi;
#if (CONICBUNDLE_DEBUG>=1)
 bool hb=A.is_init;A.is_init=B.is_init;B.is_init=hb;
#endif
}

inline Indexmatrix rand(Integer rows,Integer cols,Integer lb,Integer ub,CH_Tools::GB_rand* random_generator)
{Indexmatrix A; return A.rand(rows,cols,lb,ub,random_generator);}
inline Indexmatrix sign(const Indexmatrix& A)
          {Indexmatrix B(A); return B.sign();}

inline Indexmatrix operator>(const Indexmatrix &A,const Indexmatrix &B)
{return B<A;}
inline Indexmatrix operator>=(const Indexmatrix &A,const Indexmatrix &B)
{return B<=A;}
inline Indexmatrix operator<(Integer d,const Indexmatrix &A)
{return A>d;}
inline Indexmatrix operator>(Integer d,const Indexmatrix &A)
{return A<d;}
inline Indexmatrix operator<=(Integer d,const Indexmatrix &A)
{return A>=d;}
inline Indexmatrix operator>=(Integer d,const Indexmatrix &A)
{return A<=d;}
inline Indexmatrix operator==(Integer d,const Indexmatrix &A)
{return A==d;}
inline Indexmatrix operator!=(Integer d,const Indexmatrix &A)
{return A!=d;}

inline Indexmatrix sortindex(const Indexmatrix& vec)
{Indexmatrix ind;sortindex(vec,ind);return ind;}

inline Indexmatrix find(const Indexmatrix& A)
          {return A.find();}
inline Indexmatrix find_number(const Indexmatrix& A,Integer num)
          {return A.find_number(num);}

inline std::vector<int>& assign(std::vector<int>& vec,const Indexmatrix& A)
{
  chk_init(A);
  vec.resize(A.dim());
  for(Integer i=0;i<A.dim();i++) vec[i]=int(A(i));
  return vec;
}

inline std::vector<long>& assign(std::vector<long>& vec,const Indexmatrix& A)
{
  chk_init(A);
  vec.resize(A.dim());
  for(Integer i=0;i<A.dim();i++) vec[i]=long(A(i));
  return vec;
}

}


#ifndef CH_MATRIX_CLASSES__MATRIX_HXX
#include "matrix.hxx"
#endif


#endif
