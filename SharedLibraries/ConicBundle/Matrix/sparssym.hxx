/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/sparssym.hxx

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



#ifndef CH_MATRIX_CLASSES__SPARSSYM_HXX
#define CH_MATRIX_CLASSES__SPARSSYM_HXX

/**  @file sparssym.hxx
    @brief Header declaring the class CH_Matrix_Classes::Sparsesym for sparse symmetric matrices with Real elements
    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg

*/

#ifndef CH_MATRIX_CLASSES__SPARSMAT_HXX
#include "sparsmat.hxx"
#endif

namespace CH_Matrix_Classes {

// **************************************************************************
//                            class definition
// **************************************************************************

/**@defgroup Sparsesymgroup Sparsesym (sparse, real, symmetric, n by n)
*/
  //@{

  /** @brief %Matrix class of symmetric matrices with real values of type #Real

      Any matrix element can be indexed by (i,j) or directly by the one dimensional 
      index (i+j*nr). The latter view directly corresponds to the vec() operator 
      often used in the linear algebra literature, i.e., the matrix is 
      transformed to a vector by stacking the columns on top of each other.

      Internally the matrices colinfo, colindex, colval store the diagonal 
      and the remaining lower triangle in columnwise order. For fast matrix
      multiplication generating sparse matrices additional information is
      stored in suppind and suppcol.

      We now explain the representation of diagonal and columnwise lower triangle.
      Thereby the diagonal is treated like a separate column with negative index.
      
      Suppose that 0<=k<=nr columns of diagonal and lower triangle of the
      symmetric nr*nr matrix contain nonzero elements, 
      Then the Indexmatrix colinfo is a k by 4 matrix with one row for
      each nonzero column and the elements of row k give 
       - colinfo(k,0):   index of k-th nonzero column in the full matrix 
             (these are sorted increasingly with k, the diagonal has negative index)
       - colinfo(k,1):   number of nonzero elements in this column(/diagonal)
       - colinfo(k,2):   index of first nonzero element of this column(/diagonal)
            into colval and colindex (explained below).
       - colinfo(k,3):   index of the column within the principal submatrix
            spanned by all nonzero entries (see suppind, suppcol below)
      
      Matrix colval is a vector with number of elements equal to the
      number of nonzeros in the lower triangle (including the diagonal). 
      It stores the values of the nonzero elements of the diagonal and
      of the lower triangle in columnwise and then rowwise order, i.e., 
      nonzero element (i1,j1) with i1>=j1 appears before a different nonzero 
      element (i2,j2) with i2>=j2 iff 
      (((i1==j1)&&(i2!=j2))||(j1<j2)||((j1==j2)&&(i1<i2)).

      Indexmatrix colindex is of the same size as colval and stores for
      diagonal elemenets the row indices and for offdiagonal elements
      the rowindces minus the column index in the corresponding position. 
      The shift of the offdiagonal indices helps to speed up some 
      computations.

      Note that in multplications with a (full or sparse) matrix it is 
      beneficial to know the structure of principal submatrix that is 
      spanned by the nonzero entries, i.e., to have available the 
      indices of all columns/rows that contain nonzeros. This may differ 
      considerably from the number of columns stored in colinfo, since 
      one column of the lower triangle may cover many rows whose 
      corresponding columns have no further nonzeros in the lower triangle. 
      This information is stored in suppcol and suppind.

      The Indexmatrix suppcol is a vector of dimension equal to the
      number of nonzero columns and contains, in increasing order,
      the indices of these columns.

      The Indexmatrix suppind is arranged corresponding to colval
      and colindex, but in contrast to colindex it gives the
      row indices with respect to the principal support submatrix
      as given by suppcol.   
      
   */

class Sparsesym: protected Memarrayuser
{
    friend class Indexmatrix;
    friend class Matrix;
    friend class Symmatrix;
    friend class Sparsemat;
    
private:

    Mtype mtype;   ///< used for MatrixError templates (runtime type information was not yet existing)
    Integer nr;    ///< number rows = number columns

    //column by column representation of lower triangle
    Indexmatrix colinfo;  ///< k by 4 matrix, for nonzero columns: index (<0 for diagonal), # nonzeros, first index in colindex/colval, index in suppport submatrix
    Indexmatrix colindex; 
    Matrix colval;        ///< gives the rowindex of the element at position i, (sorted increasingly per column)
    Indexmatrix suppind;  ///< index of an element with respect to the principal submatrix spanned by the entire support
    Indexmatrix suppcol;  ///< the index of the support column in the original matrix
        
    Real tol;          ///< >0, if abs(value)<tol, then value is taken to be zero

    bool is_init;      ///< flag whether memory is initialized, it is only used if CONICBUNDLE_DEBUG is defined

    /// initialize the matrix to a 0x0 matrix without storage
    inline void init_to_zero();            

    /// removes zeros and updates suppind and suppcol
    void update_support();     

public:

  //----------------------------------------------------
  //----  constructors, destructor, and initialization
  //----------------------------------------------------

  /** @name Constructors, Destructor, and Initialization (Members)
   */
  //@{

  /// empty matrix
  inline Sparsesym(); 
  /// copy constructor, *this=d*A
  inline Sparsesym(const Sparsesym& A,Real d=1.); 
  /// initialize to zero-matrix of size nr*nr     
  inline Sparsesym(Integer nr);
  /// initialize to size nr*nr and nz nonzeros so that this(ini[i],inj[i])=val[i] for i=0,..,nz-1; specify only one of (i,j) and (j,i), multiple elements are summed up.
  inline Sparsesym(Integer nr,Integer nz,
		   const Integer *ini,const Integer *inj,const Real* va);
  /// initialize to size nr*nr and nz nonzeros so that this(ini(i),inj(i))=val[i] for i=0,..,nz-1; specify only one of (i,j) and (j,i), multiple elements are summed up.
  inline Sparsesym(Integer nr,Integer nz,
		   const Indexmatrix& ini,const Indexmatrix& inj, const Matrix& va);

#if (CONICBUNDLE_DEBUG>=1)
   /// after external initialization, call matrix.set_init(true) (not needed if CONICBUNDLE_DEBUG is undefined)
 void set_init(bool i){is_init=i;}
  /// returns true if the matrix has been declared initialized (not needed if CONICBUNDLE_DEBUG is undefined)
  int get_init() const {return is_init;} 
#else 
  /// after external initialization, call matrix.set_init(true) (not needed if CONICBUNDLE_DEBUG is undefined) 
  void set_init(bool /* i */){}
  /// returns true if the matrix has been declared initialized (not needed if CONICBUNDLE_DEBUG is undefined)
  bool get_init() const {return true;}
#endif 

  /// initialize to *this=A*d
  inline Sparsesym& init(const Sparsesym&,Real d=1.);   
  /// initialize to *this=d*(A+transpose(A))/2., abs(values)<tol are removed from the support
  inline Sparsesym& init(const Matrix&,Real d=1.);       
  /// initialize to *this=d*(A+transpose(A))/2., zeros are removed from the support
  inline Sparsesym& init(const Indexmatrix&,Real d=1.);  
  /// initialize to *this=A*d, abs(values)<tol are removed from the support
  inline Sparsesym& init(const Symmatrix&,Real d=1.); 
  /// initialize to *this=d*(A+transpose(A))/2.
  inline Sparsesym& init(const Sparsemat&,Real d=1.);    
  /// initialize to zero-matrix of size nr*nr     
  inline Sparsesym& init(Integer nr);     
  /// initialize to size nr*nr and nz nonzeros so that this(ini[i],inj[i])=val[i] for i=0,..,nz-1; specify only one of (i,j) and (j,i), multiple elements are summed up.
  Sparsesym& init(Integer nr,Integer nz,
		  const Integer *ini,const Integer *inj,const Real* va);
  /// initialize to size nr*nr and nz nonzeros so that this(ini(i),inj(i))=val[i] for i=0,..,nz-1; specify only one of (i,j) and (j,i), multiple elements are summed up.
  Sparsesym& init(Integer nr,Integer nz,
		  const Indexmatrix& ini,const Indexmatrix& inj, const Matrix& va);
  
  //initialize to the same support as A but with constant value d; the same support will be generated even for d=0. 
  Sparsesym& init_support(const Sparsesym& A,Real d=0.);

  /// set tolerance for recognizing zero values to t
  void set_tol(Real t){tol=t;}


  //@}

  /** @name Conversions from other Matrix Classes (Members)
   */
  //@{

  /// initialize to *this=d*(A+transpose(A))/2., abs(values)<tol are removed from the support
  inline Sparsesym(const Matrix&,Real d=1.);         
  /// initialize to *this=d*(A+transpose(A))/2., zeros are removed from the support
  inline Sparsesym(const Indexmatrix&,Real d=1.);    
  /// initialize to *this=A*d, abs(values)<tol are removed from the support
  inline Sparsesym(const Symmatrix&,Real d=1.);       
  /// initialize to *this=d*(A+transpose(A))/2.
  inline Sparsesym(const Sparsemat&,Real d=1.);

  //@}


  //----------------------------------------------------
  //----  size and type information
  //----------------------------------------------------

  /** @name Size and Type Information (Members)
   */
  //@{

  /// returns the number of rows in _nr and the number of columns in _nc
  void dim(Integer& r, Integer& c) const {r=nr; c=nr;}

  /// returns the dimension rows * columns when the matrix is regarded as a vector
  Integer dim() const {return nr*nr;}

  /// returns the row dimension
  Integer rowdim() const {return nr;}

  /// returns the column dimension
  Integer coldim() const {return nr;}

  /// returns the number of nonzeros in the lower triangle (including diagonal)
  Integer nonzeros() const {return colval.dim();}
  
  /// returns the type of the matrix, MTmatrix
  Mtype get_mtype() const {return mtype;}
  
  //@}

    
  //--------------------------------
  //----  Indexing and Submatrices
  //--------------------------------

  /** @name Indexing and Submatrices (Members)
   */
  //@{
    
  /// returns value of element (i,j) of the matrix (rowindex i, columnindex j)
  Real operator()(Integer i,Integer j) const;
  /// returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
  Real operator()(Integer i) const;    //vector of stacked columns

  /// returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
   Real operator[](Integer i) const;    //vector of stacked columns

  /// returns information on nozero diagonal/columns, k by 4, listing: index (<0 for diagonal), # nonzeros, first index in colindex/colval, index in suppport submatrix
  const Indexmatrix& get_colinfo() const {return colinfo;}

  /// returns the index vector of the column representation holding the row index for each element
  const Indexmatrix& get_colindex() const {return colindex;}
  /// returns the value vector of the column representation holding the value for each element
  const Matrix& get_colval() const {return colval;}
  /// returns the index vector of the column representation holding the row index w.r.t. the principal support submatrix for each element
  const Indexmatrix& get_suppind() const {return suppind;}
  /// returns the vector listing in ascending order the original column indices of the principal support submatrix
  const Indexmatrix& get_suppcol() const {return suppcol;}

  /// stores the nz nonzero values of the lower triangle of *this in I,J,val so that this(I(i),J(i))=val(i) for i=0,...,nz-1 and dim(I)=dim(J)=dim(val)=nz (ordered as in row representation)
  void get_edge_rep(Indexmatrix& I, Indexmatrix& J, Matrix& val) const;

  /// returns 1 if A is of the same dimension and the support of A is contained in the support of *this, 0 otherwise 
  int contains_support(const Sparsesym& A) const;

  /// returns 0 if (i,j) is not in the support, 1 otherwise
  int check_support(Integer i,Integer j) const;
      

  //@}


  /** @name Indexing and Submatrices (Friends)
   */
  //@{

    

  /// returns the diagonal of A as a dense Matrix vector
  friend Matrix diag(const Sparsesym& A);      //=(A(1,1),A(2,2),...)^t
  /// forms a sparse symmetrix matrix having vector A on its diagonal
  friend Sparsesym sparseDiag(const Matrix& A,Real tol=SPARSE_ZERO_TOL);

  /// swap the content of the two sparse matrices A and B (involves no copying)
  friend void swap(Sparsesym& A, Sparsesym& B);

  //@}
  

  //------------------------------
  //----  BLAS-like Routines
  //------------------------------

  /** @name BLAS-like Routines (Members)
   */
  //@{

  ///sets *this=d*A and returns *this
  Sparsesym& xeya(const Sparsesym& A,Real d=1.);
  ///sets and returns *this=d*(A+transpose(A))/2. where abs(values)<tol are removed from the support
  Sparsesym& xeya(const Matrix& A,Real d=1.);
  ///sets and returns *this=d*(A+transpose(A))/2. where zeros are removed from the support
  Sparsesym& xeya(const Indexmatrix& A,Real d=1.);

  //@}

  /** @name BLAS-like Routines (Friends)
   */
  //@{

  ///returns x= alpha*y+beta*x, where y may be transposed (ytrans=1); if beta==0. then x is initialized to the correct size
  inline friend Sparsesym& xbpeya(Sparsesym& x,const Sparsesym& y,Real alpha=1.,Real beta=0.);
  
  ///returns x= alpha*y+beta*z; x is initialized to the correct size
  friend Sparsesym& xeyapzb(Sparsesym& x,const Sparsesym& y,const Sparsesym& z,Real alpha=1.,Real beta=1.);
  
  ///returns C=beta*C+alpha*AA^T (or A^TA), but only on the current support of C
  friend Sparsesym& support_rankadd(const Matrix& A,Sparsesym& C,
				    Real alpha=1.,Real beta=0.,int trans=0);

  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to B; if beta==0. then C is initialized to the correct size
  friend Matrix& genmult(const Sparsesym& A,const Matrix& B,Matrix &C,
			 Real alpha=1.,Real beta=0., int btrans=0);

  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to B; if beta==0. then C is initialized to the correct size
  friend Matrix& genmult(const Matrix& A,const Sparsesym& B,Matrix &C,
			 Real alpha=1.,Real beta=0., int atrans=0);

  //@}
  

  //------------------------------
  //----  usual operators
  //------------------------------

  /** @name Usual Arithmetic Operators (Members)
   */
  //@{

  ///
  inline Sparsesym& operator=(const Sparsesym &A);
  ///
  inline Sparsesym& operator+=(const Sparsesym &v);
  ///
  inline Sparsesym& operator-=(const Sparsesym &v);
  ///
  inline Sparsesym operator-() const;
 
  ///
  inline Sparsesym& operator*=(Real d);
  /// ATTENTION: d is NOT checked for 0
  inline Sparsesym& operator/=(Real d);
  
  ///sets *this=(A+transpose(A))/2. removing abs(values)<tol; returns *this
  inline Sparsesym& operator=(const Matrix &A);
  ///sets *this=(A+transpose(A))/2. removing zeros; returns *this
  inline Sparsesym& operator=(const Indexmatrix &A);

  ///transposes itself (at almost no cost)  
  Sparsesym& transpose(){return *this;}            

   
  //@}

  /** @name Usual Arithmetic Operators (Friends)
   */
  //@{
  
  ///
  friend inline Sparsesym operator+(const Sparsesym &A,const Sparsesym& B);
  ///
  friend inline Sparsesym operator-(const Sparsesym &A,const Sparsesym& B);
  ///
  friend inline Sparsesym operator*(const Sparsesym& A,Real d);
  ///
  friend inline Sparsesym operator*(Real d,const Sparsesym &A);
  /// ATTENTION: d is NOT checked for 0
  friend inline Sparsesym operator/(const Sparsesym& A,Real d);
  
  
  ///
  friend inline Matrix operator*(const Matrix& A,const Sparsesym& B);
  ///
  friend inline Matrix operator*(const Sparsesym& A,const Matrix& B);
  ///
  friend inline Matrix operator+(const Matrix& A,const Sparsesym& B);
  ///
  friend inline Matrix operator+(const Sparsesym& A,const Matrix& B);
  ///
  friend inline Matrix operator-(const Matrix& A,const Sparsesym& B);
  ///
  friend inline Matrix operator-(const Sparsesym& A,const Matrix& B);
  
  
  /// (drop it or use a constructor instead)
  friend inline Sparsesym transpose(const Sparsesym& A);

  //@}

  //------------------------------------------
  //----  Connections to other Matrix Classes
  //------------------------------------------
  
  /** @name Connections to other Classes (Members)
   */
  //@{

  /// sets and returns *this=A*d where abs(values)<tol are removed from the support
  Sparsesym& xeya(const Symmatrix& A,Real d=1.);
  /// sets and returns *this=d.*(A+transpose(A))/2.
  Sparsesym& xeya(const Sparsemat& A,Real d=1.);

  /// sets and returns *this=A where abs(values)<tol are removed from the support
  inline Sparsesym& operator=(const Symmatrix &A);
  ///
  inline Symmatrix operator+(const Symmatrix &A) const;
  ///
  inline Symmatrix operator-(const Symmatrix &A) const;   

  /// sets and returns *this=(A+transpose(A))/2.
  inline Sparsesym& operator=(const Sparsemat &A);
  /// compute (*this)*A and return the result in a Sparsemat
  Sparsemat sparsemult(const Matrix& A) const;

  //@}

  /** @name Connections to other Classes (Friends)
   */
  //@{

  ///
  friend inline Symmatrix operator+(const Sparsesym &A,const Symmatrix &B);
  ///
  friend inline Symmatrix operator+(const Symmatrix &A,const Sparsesym &B);
  ///
  friend inline Symmatrix operator-(const Sparsesym& A,const Symmatrix &B);
  ///
  friend inline Symmatrix operator-(const Symmatrix &A,const Sparsesym &B);
  ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  friend Real ip(const Symmatrix& A, const Sparsesym& B); //=trace(B^t*A)
  ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  friend inline Real ip(const Sparsesym& A, const Symmatrix& B);

  ///returns C=beta*C+alpha*A*B, where B may be transposed; if beta==0. then C is initialized to the correct size
  friend Matrix& genmult(const Sparsesym& A,const Sparsemat& B,Matrix &C,
			 Real alpha=1.,Real beta=0.,int btrans=0);
  ///returns C=beta*C+alpha*A*B, where A may be transposed; if beta==0. then C is initialized to the correct size
  friend Matrix& genmult(const Sparsemat& A,const Sparsesym& B,Matrix &C,
			 Real alpha=1.,Real beta=0.,int atrans=0);
  ///
  friend inline Matrix operator*(const Sparsesym& A,const Sparsemat& B);
  ///
  friend inline Matrix operator*(const Sparsemat& A,const Sparsesym& B);

  //@}

  //------------------------------
  //----  Elementwise Operations
  //------------------------------

  /** @name Elementwise Operations (Friends)
   */
  //@{

  /// sets (*this)(i,j)=abs((*this)(i,j)) for all i,j and returns *this
  friend Sparsesym abs(const Sparsesym& A);                 

  //@}

  //----------------------------
  //----  Numerical Methods
  //----------------------------

  /** @name Numerical Methods (Friends)
   */
  //@{

  /// returns the sum of the diagonal elements A(i,i) over all i
  friend Real trace(const Sparsesym& A);               //=sum(diag(A))
  ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  friend Real ip(const Sparsesym& A, const Sparsesym& B); //=trace(B^t*A)
  ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  friend Real ip(const Matrix& A, const Sparsesym& B); //=trace(B^t*A)
  ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  friend inline Real ip(const Sparsesym& A, const Matrix& B); 
  ///returns the Frobenius norm of A, i.e., the square root of the sum of A(i,j)*A(i,j) over all i,j 
  friend Real norm2(const Sparsesym& A); //=sqrt(ip(A,A))

  ///returns a row vector holding the sum over all rows, i.e., (1 1 ... 1)*A
  friend Matrix sumrows(const Sparsesym& A);   //=(1 1 1 ... 1)*A
  ///returns a column vector holding the sum over all columns, i.e., A*(1 1 ... 1)^T
  friend inline Matrix sumcols(const Sparsesym& A);    //=A*(1 1 1 ... 1)^t
  ///returns the sum over all elements of A, i.e., (1 1 ... 1)*A*(1 1 ... 1)^T
  friend Real sum(const Sparsesym& A);         //=(1 1 ... 1)*A*(1 1 ... 1)^t

  //@}


  //---------------------------------------------
  //----  Comparisons / Max / Min / sort / find
  //---------------------------------------------

  /** @name Equal (Members)
   */
  //@{

  /// returns 1 if both matrices are identical, 0 otherwise
   friend int equal(const Sparsesym& A, const Sparsesym& B,Real eqtol);

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
		 int precision=0,   ///< number of most significant digits, default=4
		 int width=0,       ///< field width, default = precision+6
		 int screenwidth=0  ///< maximum number of characters in one output line, default = 80
		 ) const;

  //@}

  /** @name Input, Output (Friends)
   */
  //@{

  ///output format (lower triangle): nr nz \\n i1 j1 val1\\n i2 j2 val2\\n ... inz jnz valnz\\n 
  friend std::ostream& operator<<(std::ostream& o,const Sparsesym &v);
  ///input format (lower triangle): nr nz \\n i1 j1 val1\\n i2 j2 val2\\n ... inz jnz valnz\\n 
  friend std::istream& operator>>(std::istream& i,Sparsesym &v);

  //@}

};

//@}

// **************************************************************************
//                make non inline friends available outside
// **************************************************************************

  Matrix diag(const Sparsesym& A);      //=(A(1,1),A(2,2),...)^t
  Sparsesym sparseDiag(const Matrix& A,Real tol);
  void swap(Sparsesym& A, Sparsesym& B);
  Sparsesym& xeyapzb(Sparsesym& x,const Sparsesym& y,
		     const Sparsesym& z,Real alpha,Real beta);
  Sparsesym& support_rankadd(const Matrix& A,Sparsesym& C,
			     Real alpha,Real beta,int trans);
  Matrix& genmult(const Sparsesym& A,const Matrix& B,Matrix &C,
		  Real alpha,Real beta, int btrans);
  Matrix& genmult(const Matrix& A,const Sparsesym& B,Matrix &C,
		  Real alpha,Real beta, int atrans);
  Real ip(const Symmatrix& A, const Sparsesym& B); //=trace(B^t*A)
  Matrix& genmult(const Sparsesym& A,const Sparsemat& B,Matrix &C,
		  Real alpha,Real beta,int btrans);
  Matrix& genmult(const Sparsemat& A,const Sparsesym& B,Matrix &C,
		  Real alpha,Real beta,int atrans);
  Sparsesym abs(const Sparsesym& A);
  Real trace(const Sparsesym& A);               //=sum(diag(A))
  Real ip(const Sparsesym& A, const Sparsesym& B); //=trace(B^t*A)
  Real ip(const Matrix& A, const Sparsesym& B); //=trace(B^t*A)
  Real norm2(const Sparsesym& A); //=sqrt(ip(A,A))
  Matrix sumrows(const Sparsesym& A);   //=(1 1 1 ... 1)*A
  Real sum(const Sparsesym& A);         //=(1 1 ... 1)*A*(1 1 ... 1)^t
  std::ostream& operator<<(std::ostream& o,const Sparsesym &v);
  std::istream& operator>>(std::istream& i,Sparsesym &v);

// **************************************************************************
//                   implementation of inline functions
// **************************************************************************

inline void Sparsesym::init_to_zero()
{
 mtype=MTsparse;
 nr=0;
 tol=SPARSE_ZERO_TOL;
#if (CONICBUNDLE_DEBUG>=1)
 is_init=1;
#endif
}

inline Sparsesym& Sparsesym::init(const Sparsesym& A,Real d)
{ return xeya(A,d); }

inline Sparsesym& Sparsesym::init(const Matrix& A,Real d)
{ return xeya(A,d);}

inline Sparsesym& Sparsesym::init(const Indexmatrix& A,Real d)
{ return xeya(A,d);}

inline Sparsesym& Sparsesym::init(Integer r)
{
nr=r;
colinfo.init(0,0,Integer(0)); colindex.init(0,0,Integer(0)); colval.init(0,0,0.);
chk_set_init(*this,1);
return *this;
}

inline Sparsesym::Sparsesym()
{ init_to_zero(); chk_set_init(*this,1);}
inline Sparsesym::Sparsesym(Integer r)
{ init_to_zero(); init(r);}
inline Sparsesym::Sparsesym(const Sparsesym& A,Real d):Memarrayuser()
{ init_to_zero(); xeya(A,d);}
inline Sparsesym::Sparsesym(const Matrix& A,Real d)
{ init_to_zero(); xeya(A,d);}
inline Sparsesym::Sparsesym(const Indexmatrix& A,Real d)
{ init_to_zero(); xeya(A,d);}

inline Sparsesym::Sparsesym(Integer in_nr,Integer nz,
       const Integer *ini,const Integer *inj,const Real* va)
{ init_to_zero(); init(in_nr,nz,ini,inj,va);}

inline Sparsesym::Sparsesym(Integer in_nr,Integer nz,
              const Indexmatrix& ini,const Indexmatrix& inj, const Matrix& va)
{ init_to_zero(); init(in_nr,nz,ini,inj,va);}

inline Real Sparsesym::operator()(Integer i) const
{ return (*this)(i%nr,i/nr); }
inline Real Sparsesym::operator[](Integer i) const
{ return (*this)(i%nr,i/nr); }

inline Sparsesym& xbpeya(Sparsesym& x,const Sparsesym& y,Real alpha,Real beta)
{
  if (beta==0.) return x.init(y,alpha); 
  Sparsesym B; xeyapzb(B,x,y,beta,alpha); swap(x,B); return x;
}

inline Sparsesym& Sparsesym::operator=(const Sparsesym &A)
{ return xeya(A);} 
inline Sparsesym& Sparsesym::operator+=(const Sparsesym &A)
{ Sparsesym B; xeyapzb(B,*this,A); swap(*this,B); return *this;}
inline Sparsesym& Sparsesym::operator-=(const Sparsesym &A)
{ Sparsesym B; xeyapzb(B,*this,A,1.,-1.); swap(*this,B); return *this;}
inline Sparsesym& Sparsesym::operator*=(Real d)
{chk_init(*this); colval*=d; return *this;}
inline Sparsesym& Sparsesym::operator/=(Real d)
{chk_init(*this);colval/=d;return *this;}
inline Sparsesym Sparsesym::operator-() const
{ return Sparsesym(*this,-1.); }

    
inline Sparsesym& Sparsesym::operator=(const Matrix &A)
{ return xeya(A);}
inline Sparsesym& Sparsesym::operator=(const Indexmatrix &A)
{ return xeya(A);}

inline Sparsesym operator+(const Sparsesym& A,const Sparsesym& B) 
{Sparsesym C; return xeyapzb(C,A,B);}
inline Sparsesym operator-(const Sparsesym& A,const Sparsesym& B) 
{Sparsesym C; return xeyapzb(C,A,B,1.,-1.);}

inline Sparsesym operator*(const Sparsesym &A,Real d) 
    { return Sparsesym(A,d);}
inline Sparsesym operator*(Real d,const Sparsesym &A)
    { return Sparsesym(A,d);}
inline Sparsesym operator/(const Sparsesym &A,Real d)
    { return Sparsesym(A,1/d);}

inline Matrix operator*(const Sparsesym& A,const Matrix& B)
{ Matrix C; return genmult(A,B,C); }
inline Matrix operator*(const Matrix& A,const Sparsesym& B)
{ Matrix C; return genmult(A,B,C); }
inline Matrix operator+(const Matrix& A,const Sparsesym& B)
{ Matrix C(A); return C.xpeya(B);}
inline Matrix operator+(const Sparsesym& A,const Matrix& B)
{ Matrix C(B); return C.xpeya(A);}
inline Matrix operator-(const Matrix& A,const Sparsesym& B)
{ Matrix C(A); return C.xpeya(B,-1.);}
inline Matrix operator-(const Sparsesym& A,const Matrix& B)
{ Matrix C(B,-1.); return C.xpeya(A);}

inline Real ip(const Sparsesym& A, const Matrix& B)
{return ip(B,A);}
inline Sparsesym transpose(const Sparsesym& A)
{return Sparsesym(A);}

inline Matrix sumcols(const Sparsesym& A)    //=A*(1 1 1 ... 1)^t
{Matrix s(sumrows(A)); return s.transpose();}


inline Matrix::Matrix(const Sparsesym& A,Real d)

{init_to_zero(); xeya(A,d); }

inline Matrix& Matrix::init(const Sparsesym& A,Real d)

{return xeya(A,d); }

inline Matrix& Matrix::operator=(const Sparsesym& A)

{return xeya(A); }

inline Matrix& Matrix::operator*=(const Sparsesym& A)

{ Matrix C; return *this=genmult(*this,A,C);}

inline Matrix& Matrix::operator+=(const Sparsesym& A)

{return xpeya(A); }

inline Matrix& Matrix::operator-=(const Sparsesym& A)

{return xpeya(A,-1.); }


inline Sparsesym::Sparsesym(const Symmatrix& A,Real d)
{ init_to_zero(); xeya(A,d);}
inline Sparsesym& Sparsesym::init(const Symmatrix& A,Real d)
{ return xeya(A,d);}
inline Sparsesym& Sparsesym::operator=(const Symmatrix& A)
{ return xeya(A);}
inline Symmatrix Sparsesym::operator+(const Symmatrix &A) const
{Symmatrix B(A); B+=*this; return B;}
inline Symmatrix Sparsesym::operator-(const Symmatrix &A) const
{Symmatrix B(A); B-=*this; return B;}
inline Symmatrix operator+(const Sparsesym &A,const Symmatrix &B)
{Symmatrix C(B);return C.xpeya(A);}
inline Symmatrix operator+(const Symmatrix &A,const Sparsesym &B)
{Symmatrix C(A);return C.xpeya(B);}
inline Symmatrix operator-(const Sparsesym& A,const Symmatrix &B)
{Symmatrix C(B,-1.);return C.xpeya(A);}
inline Symmatrix operator-(const Symmatrix &A,const Sparsesym &B)
{Symmatrix C(A);return C.xpeya(B,-1);}
inline Real ip(const Sparsesym& A, const Symmatrix& B)
{ return ip(B,A);}    
inline Sparsesym::Sparsesym(const Sparsemat& A,Real d)
    { init_to_zero(); xeya(A,d);}
inline Sparsesym& Sparsesym::init(const Sparsemat& A,Real d)
    { return xeya(A,d);}
inline Sparsesym& Sparsesym::operator=(const Sparsemat& A)
    { return xeya(A);}
inline Matrix operator*(const Sparsesym& A,const Sparsemat& B)
    { Matrix C; return genmult(A,B,C);} 
inline Matrix operator*(const Sparsemat& A,const Sparsesym& B)
    { Matrix C; return genmult(A,B,C);} 

inline Sparsemat::Sparsemat(const Sparsesym& A,Real d)
{ init_to_zero(); xeya(A,d);}
inline Sparsemat& Sparsemat::init(const Sparsesym& A,Real d) 
{ return xeya(A,d);}   
inline Sparsemat& Sparsemat::operator=(const Sparsesym &A)
{ return xeya(A);}  
inline Matrix operator*(const Symmatrix& A,const Sparsemat& B)
{ Matrix C; return genmult(A,B,C,1.,0.,0);} 
inline Matrix operator*(const Sparsemat& A,const Symmatrix& B)
{ Matrix C; return genmult(A,B,C,1.,0.,0);}  

inline Symmatrix::Symmatrix(const Sparsesym& A,Real d)
    { init_to_zero(); xeya(A,d);}
inline Symmatrix& Symmatrix::init(const Sparsesym& A,Real d)          
    { return xeya(A,d);}
inline Symmatrix& Symmatrix::operator=(const Sparsesym& A)                
    { return xeya(A); }
inline Symmatrix& Symmatrix::operator+=(const Sparsesym& A)               
    { return xpeya(A);}
inline Symmatrix& Symmatrix::operator-=(const Sparsesym& A)           
    { return xpeya(A,-1.);}

/** @name Equal (Members)
 */
//@{

/// returns 1 if both matrices are identical, 0 otherwise
int equal(const Sparsesym& A, const Sparsesym& B,Real eqtol=1e-10);


//@}    

}


#endif

