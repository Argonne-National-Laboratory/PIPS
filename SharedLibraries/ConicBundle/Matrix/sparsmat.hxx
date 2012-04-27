/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/sparsmat.hxx

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



#ifndef CH_MATRIX_CLASSES__SPARSMAT_HXX
#define CH_MATRIX_CLASSES__SPARSMAT_HXX

/**  @file sparsmat.hxx
    @brief Header declaring the class CH_Matrix_Classes::Sparsemat for sparse matrices with Real elements
    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg

*/


#ifndef CH_MATRIX_CLASSES__SYMMAT_HXX
#include "symmat.hxx"
#endif

namespace CH_Matrix_Classes {

//everything involving a "Sparsesym" is implemented in "sparssym.cxx/hxx"

#define SPARSE_ZERO_TOL 1e-60
/**@defgroup Sparsematgroup Sparsemat (sparse, real, n by m)
*/
  //@{

  /** @brief %Matrix class of sparse matrices with real values of type #Real

      Any matrix element can be indexed by (i,j) or directly by the one dimensional 
      index (i+j*nr). The latter view directly corresponds to the vec() operator 
      often used in the linear algebra literature, i.e., the matrix is 
      transformed to a vector by stacking the columns on top of each other.

      Internally a sparse matrix of size nr x nc is automatically stored 
      in sparse columnwise as well as rowwise format by means of additional
      objects of type Indexmatrix and Matrix. The purpose of this somewhat
      complicated format is to allow fast multiplication from left and
      right without increasing the amount of information stored by more 
      than a constant factor times the number of nonzero elements.

      We now explain the columnwise format (the rowwise is structured
      in the same manner). It is given by the matrices colinfo, colval,
      and colindex. 
   
      Suppose that 0<=k<=nc columns of the matrix contain nonzero elements, 
      then the Indexmatrix colinfo is a k by 3 matrix with one row for
      each nonzero column and the elements of row k give 
       - colinfo(k,0):   index of k-th nonzero column in the full matrix 
             (these are sorted increasingly with k)
       - colinfo(k,1):   number of nonzero elements in this column
       - colinfo(k,2):   index of first nonzero element of this column
            into colval and colindex (explained below).

      Matrix colval is a vector with number of elements equal to the
      number of nonzeros in the entire matrix. It stores the values of 
      the nonzero elements in columnwise and then rowwise order, i.e., 
      nonzero element (i1,j1) appears before a different nonzero element 
      (i2,j2) iff ((j1<j2)||((j1==j2)&&(i1<i2)).

      Indexmatrix colindex is of the same size as colval and stores the
      corresponding row indices in the corresponding position. 

      Finding an element (i,j) is thus done as follows: First find
      the column index j in the first column of colinfo by binary search
      yielding a colinfo-rowindex k if successful; then find the row index 
      i within colindex[colinfo(k,2)]...colindex[colinfo(k,2)+colinfo(k,1)-1]
      again by binary search.  
               

   */

class Sparsemat: protected Memarrayuser
{
    friend class Indexmatrix;
    friend class Matrix;
    friend class Symmatrix;
    friend class Sparsesym;
    
private:

    Mtype mtype;   ///< used for MatrixError templates (runtime type information was not yet existing)
    Integer nr,    ///< number of rows
            nc;    ///< number of columns      

    //column by column representation
    Indexmatrix colinfo;  ///< k by 3, for nonzero columns: index, # nonzeros, first index in colindex/colval
    Indexmatrix colindex; ///< gives the rowindex of the element at position i, (sorted increasingly per column)
    Matrix colval;        ///< gives the value of the element at position i

    //row by row representation (see column by column)
    Indexmatrix rowinfo; ///< k by 3, for nonzero rows: index, # nonzeros, first index in colindex/colval
    Indexmatrix rowindex;///< gives the column index of the element at position i, (sorted increasingly per row)
    Matrix rowval;       ///< gives the value of the element at position i

        
    Real tol;          ///< >0, if abs(value)<tol, then value is taken to be zero

    bool is_init;      ///< flag whether memory is initialized, it is only used if CONICBUNDLE_DEBUG is defined

    /// initialize the matrix to a 0x0 matrix without storage
    void init_to_zero();          

public:

  //----------------------------------------------------
  //----  constructors, destructor, and initialization
  //----------------------------------------------------

  /** @name Constructors, Destructor, and Initialization (Members)
   */
  //@{

  /// empty matrix
  inline Sparsemat();                      
  /// copy constructor, *this=d*A, abs(values)<tol are removed from the support
  inline Sparsemat(const Sparsemat& A,Real d=1.);       
  /// initialize to zero-matrix of size nr*nc     
  inline Sparsemat(Integer nr,Integer nc);  
  /// initialize to size nr*nc and nz nonzeros so that this(ini[i],inj[i])=val[i] for i=0,..,nz-1; multiple elements are summed up.
  inline Sparsemat(Integer nr,Integer nc,Integer nz,
		   const Integer *ini,const Integer *inj,const Real* va);
  /// initialize to size nr*nc and nz nonzeros so that this(ini(i),inj(i))=val(i) for i=0,..,nz-1; multiple elements are summed up.
  inline Sparsemat(Integer nr,Integer nc,Integer nz,
		   const Indexmatrix& ini,const Indexmatrix& inj, const Matrix& va);

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
  inline Sparsemat& init(const Sparsemat& A,Real d=1.);       
  /// initialize to *this=A*d, abs(values)<tol are removed from the support
  inline Sparsemat& init(const Matrix& A,Real d=1.);          
  /// initialize to *this=A*d, zeros are removed from the support
  inline Sparsemat& init(const Indexmatrix& A,Real d=1.);     
  /// initialize to *this=A*d, abs(values)<tol are removed from the support
  inline Sparsemat& init(const Symmatrix&,Real d=1.);     
  /// initialize to *this=A*d
  inline Sparsemat& init(const Sparsesym&,Real d=1.);    
  /// initialize to zero-matrix of size nr*nc     
  inline Sparsemat& init(Integer nr,Integer nc);  //zero-matrix of size nr*nc
  /// initialize to size nr*nc and nz nonzeros so that this(ini[i],inj[i])=val[i] for i=0,..,nz-1; multiple elements are summed up.
  Sparsemat& init(Integer nr,Integer nc,Integer nz,
		  const Integer *ini,const Integer *inj,const Real* va);
  /// initialize to size nr*nc and nz nonzeros so that this(ini(i),inj(i))=val(i) for i=0,..,nz-1; multiple elements are summed up.
  Sparsemat& init(Integer nr,Integer nc,Integer nz,
		  const Indexmatrix& ini,const Indexmatrix& inj, const Matrix& va);

  /// set tolerance for recognizing zero values to t
  void set_tol(Real t){tol=t;}

  //@}

  /** @name Conversions from other Matrix Classes (Members)
   */
  //@{

  /// initialize to *this=d*A, abs(values)<tol are removed from the support
  inline Sparsemat(const Matrix& A, Real d=1.);          
  /// initialize to *this=d*A, zeros are removed from the support
  inline Sparsemat(const Indexmatrix& A,Real d=1.);
  /// initialize to *this=d*A, abs(values)<tol are removed from the support
  inline Sparsemat(const Symmatrix& A,Real d=1.); 
  /// initialize to *this=d*A
  inline Sparsemat(const Sparsesym& A,Real d=1.);

  //@}

  //----------------------------------------------------
  //----  size and type information
  //----------------------------------------------------

  /** @name Size and Type Information (Members)
   */
  //@{

  /// returns the number of rows in _nr and the number of columns in _nc
  void dim(Integer& in_nr, Integer& in_nc) const {in_nr=nr; in_nc=nc;}

  /// returns the dimension rows * columns when the matrix is regarded as a vector
  Integer dim() const {return nr*nc;}

  /// returns the row dimension
  Integer rowdim() const {return nr;}

  /// returns the column dimension
  Integer coldim() const {return nc;}

  /// returns the number of nonzeros
  Integer nonzeros() const {return colval.dim();}

  /// returns the number of nonzeros in column i; if nonzeros>0 and startind!=0 then the index of the first nonzero in colindex/colval is stored there
  Integer col_nonzeros(Integer i,Integer* startind=0) const;

  /// returns the number of nonzeros in row i; if nonzeros>0 and startind!=0 then the index of the first nonzero in rowindex/rowval is stored there
  Integer row_nonzeros(Integer i,Integer* startind=0) const;

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
  inline Real operator()(Integer i) const;           

  /// returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
  inline Real operator[](Integer i) const; 

   /// returns column i copied to a new sparse matrix 
  Sparsemat col(Integer i) const;

  /// returns row i copied to a new sparse matrix
  Sparsemat row(Integer i) const;
  /// returns a sparse matrix of size this->rowdim() x vec.dim(), with column i a copy of column vec(i) of *this
  Sparsemat cols(const Indexmatrix& ind) const;
  /// returns a sparse matrix of size vec.dim() x this->rowdim(), with row i a copy of row vec(i) of *this
  Sparsemat rows(const Indexmatrix& ind) const;
  
  /// all rows indexed by vector ind are deleted, no row should appear twice in ind, remaining rows are moved up keeping their order, returns *this
  Sparsemat& delete_rows(const Indexmatrix& ind);
  /// all colmuns indexed by vector ind are deleted, no column should appear twice in ind, remaining columns are moved up keeping their order, returns *this
  Sparsemat& delete_cols(const Indexmatrix& ind);

  /// insert the row vector v before row i, 0<=i<= row dimension, for i==row dimension the row is appended below; appending to a 0x0 matrix is allowed, returns *this
  Sparsemat& insert_row(Integer i,const Sparsemat& v); 
  /// insert a column before column i, 0<=i<= column dimension, for i==column dimension the column is appended at the right; appending to a 0x0 matrix is allowed, returns *this
  Sparsemat& insert_col(Integer i,const Sparsemat& v); 

  /// concats sparse matrix A to the right of *this, A or *this may be the 0x0 matrix initally, returns *this
  Sparsemat& concat_right(const Sparsemat& A);
  /// concats sparse matrix A to the bottom of *this, A or *this may be the 0x0 matrix initally, returns *this
  Sparsemat& concat_below(const Sparsemat& A);

  /// returns information on nonzero columns, k by 3, listing: index, %#nonzeros, first index in colindex/colval
  const Indexmatrix& get_colinfo() const {return colinfo;}
  /// returns the index vector of the column representation holding the row index for each element
  const Indexmatrix& get_colindex() const {return colindex;}
  /// returns the value vector of the column representation holding the value for each element
  const Matrix& get_colval() const {return colval;}

  /// returns information on nonzero rows, k by 3, listing: index, %#nonzeros, first index in rowindex/rowval
  const Indexmatrix& get_rowinfo() const {return rowinfo;}
  /// returns the index vector of the row representation holding the column index for each element
  const Indexmatrix& get_rowindex() const {return rowindex;}
  /// returns the value vector of the row representation holding the value for each element
  const Matrix& get_rowval() const {return rowval;}
  
  /// stores the nz nonzero values of *this in I,J,val so that this(I(i),J(i))=val(i) for i=0,...,nz-1 and dim(I)=dim(J)=dim(val)=nz (ordered as in row representation)
  void get_edge_rep(Indexmatrix& I, Indexmatrix& J, Matrix& val) const;
  /// stores element i of the get_edge_rep() function (ordered as in row representation); returns 1 if i is out of range, 0 otherwise.
  int get_edge(Integer i,Integer& indi,Integer& indj, Real& val) const;
  /// returns 1 if A is of the same dimension and the support of A is contained in the support of *this, 0 otherwise 
  int contains_support(const Sparsemat& A) const;


  //@}


  /** @name Indexing and Submatrices (Friends)
   */
  //@{

  /// returns a new sparse matrix [A, B], i.e., it concats matrices A and B rowwise; A or B may be a 0x0 matrix
  friend inline Sparsemat concat_right(const Sparsemat& A,const Sparsemat& B);
  /// returns a new sparse matrix [A; B], i.e., it concats matrices A and B columnwise; A or B may be a 0x0 matrix
  friend inline Sparsemat concat_below(const Sparsemat& A,const Sparsemat& B);

  /// swap the content of the two sparse matrices A and B (involves no copying)
  friend void swap(Sparsemat& A, Sparsemat& B);

  //@}

  //------------------------------
  //----  BLAS-like Routines
  //------------------------------

  /** @name BLAS-like Routines (Members)
   */
  //@{

  ///sets *this=d*A and returns *this
  Sparsemat& xeya(const Sparsemat& A,Real d=1.);
  ///sets *this=d*A removing abs(values)<tol; returns *this
  Sparsemat& xeya(const Matrix& A,Real d=1.);
  ///sets *this=d*A removing zeros; returns *this
  Sparsemat& xeya(const Indexmatrix& A,Real d=1.);

  //@}

  /** @name BLAS-like Routines (Friends)
   */
  //@{

  ///returns x= alpha*y+beta*x, where y may be transposed (ytrans=1); if beta==0. then x is initialized to the correct size
  friend Sparsemat& xbpeya(Sparsemat& x,const Sparsemat& y,Real alpha=1.,Real beta=0.,int ytrans=0);
  
  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
  friend Matrix& genmult(const Sparsemat& A,const Matrix& B,Matrix &C,
			 Real alpha=1.,Real beta=0., int atrans=0,int btrans=0);
      //C=beta*C+alpha A B

  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
  friend Matrix& genmult(const Matrix& A,const Sparsemat& B,Matrix &C,
			 Real alpha=1.,Real beta=0., int atrans=0,int btrans=0);
      //C=beta*C+alpha A B
  
  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
  friend Matrix& genmult(const Sparsemat& A,const Sparsemat& B,Matrix &C,
			 Real alpha=1.,Real beta=0., int atrans=0,int btrans=0);
      //C=beta*C+alpha A B

  //@}

  //------------------------------
  //----  usual operators
  //------------------------------

  /** @name Usual Arithmetic Operators (Members)
   */
  //@{

  ///
  inline Sparsemat& operator=(const Sparsemat& A);
  ///
  inline Sparsemat& operator+=(const Sparsemat& A);
  ///
  inline Sparsemat& operator-=(const Sparsemat& A);
  ///
  inline Sparsemat operator-() const;
  
  ///
  inline Sparsemat& operator*=(Real d);
  ///
  inline Sparsemat& operator/=(Real d);
  
  ///sets *this=A removing abs(values)<tol; returns *this
  Sparsemat& operator=(const Matrix& A);
  ///sets *this=A removing zeros; returns *this
  Sparsemat& operator=(const Indexmatrix& A);

  ///transposes itself (swaps row and column representations, thus cheap)
  Sparsemat& transpose();          
   
  //@}

  /** @name Usual Arithmetic Operators (Friends)
   */
  //@{
  
  ///
  friend Sparsemat operator*(const Sparsemat& A,const Sparsemat& B);
  ///
  friend inline Sparsemat operator+(const Sparsemat &A,const Sparsemat& B);
  ///
  friend inline Sparsemat operator-(const Sparsemat &A,const Sparsemat& B);

  ///
  friend inline Sparsemat operator*(const Sparsemat& A,Real d);
  ///
  friend inline Sparsemat operator*(Real d,const Sparsemat& A);
  /// ATTENTION: d is NOT checked for 0
  friend inline Sparsemat operator/(const Sparsemat& A,Real d);
  
  ///
  friend inline Matrix operator*(const Sparsemat& A,const Matrix& B);
  ///
  friend inline Matrix operator*(const Matrix& A,const Sparsemat& B);
  ///
  friend inline Matrix operator+(const Sparsemat& A,const Matrix& B);
  ///
  friend inline Matrix operator+(const Matrix& A,const Sparsemat& B);
  ///
  friend inline Matrix operator-(const Sparsemat& A,const Matrix& B);
  ///
  friend inline Matrix operator-(const Matrix& A,const Sparsemat& B);

  ///
  friend inline Sparsemat transpose(const Sparsemat& A);

  //@}

  //------------------------------------------
  //----  Connections to other Matrix Classes
  //------------------------------------------
  
  /** @name Connections to other Classes (Members)
   */
  //@{

  /// sets *this=A*d, abs(values)<tol are removed from the support, and returns *this,
  Sparsemat& xeya(const Sparsesym& A,Real d=1.);
  /// sets *this=A, abs(values)<tol are removed from the support, returns *this
  inline Sparsemat& operator=(const Symmatrix& A);
  /// sets *this=A and returns *this
  inline Sparsemat& operator=(const Sparsesym& A);

  //@}

  /** @name Connections to other Classes (Friends)
   */
  //@{

  ///returns C=beta*C+alpha*A*B, where B may be transposed; if beta==0. then C is initialized to the correct size
  friend Matrix& genmult(const Symmatrix& A,const Sparsemat& B,Matrix& C,
			 Real alpha,Real beta,int btrans);
      //returns C=beta*C+alpha*A*B, where A and B may be transposed

  ///returns C=beta*C+alpha*A*B, where A may be transposed; if beta==0. then C is initialized to the correct size
  friend Matrix& genmult(const Sparsemat& A,const Symmatrix& B,Matrix& C,
			 Real alpha,Real beta, int atrans);

  ///returns C=beta*C+alpha*A*B, where B may be transposed; if beta==0. then C is initialized to the correct size
  friend Matrix& genmult(const Sparsesym& A,const Sparsemat& B,Matrix &C,
			 Real alpha,Real beta, int btrans);

  ///returns C=beta*C+alpha*A*B, where A may be transposed; if beta==0. then C is initialized to the correct size
  friend Matrix& genmult(const Sparsemat& A,const Sparsesym& B,Matrix &C,
			 Real alpha,Real beta,int atrans);

  /// returns C=beta*C+alpha* A*A^T, where A may be transposed; if beta==0. then C is initialized to the correct size
  friend Symmatrix& rankadd(const Sparsemat& A,Symmatrix& C,
			    Real alpha,Real beta,int trans);
      // returns C=beta*C+alpha* A*A^T, where A may be transposed
    
  /// returns C=beta*C+alpha*(A*B^T+B*A^T)/2 [or for transposed (A^T*B+B^T*A)/2]. If beta==0. then C is initiliazed to the correct size.
  friend Symmatrix& rank2add(const Sparsemat& A,const Matrix& B,Symmatrix& C,
			     Real alpha,Real beta,int trans);
      // returns C=beta*C+alpha* sym(A*B^T) [or sym(A^T*B)]

  ///
  friend inline Matrix operator*(const Symmatrix& A,const Sparsemat& B);
  ///
  friend inline Matrix operator*(const Sparsemat& A,const Symmatrix& B);

  //@}

  //------------------------------
  //----  Elementwise Operations
  //------------------------------

  /** @name Elementwise Operations (Friends)
   */
  //@{

  /// sets (*this)(i,j)=abs((*this)(i,j)) for all i,j and returns *this
  friend Sparsemat abs(const Sparsemat& A);                 

  //@}

  //----------------------------
  //----  Numerical Methods
  //----------------------------

  /** @name Numerical Methods (Members)
   */
  //@{

  ///scales each row i of (*this) by vec(i), i.e., (*this)=diag(vec)*(*this), and returns (*this)
  Sparsemat& scale_rows(const Matrix& vec);    //A=diag(vec)*A
  ///scales each column i of (*this) by vec(i), i.e., (*this)=(*this)*diag(vec), and returns (*this)
  Sparsemat& scale_cols(const Matrix& vec);    //A=A*diag(vec)


  //@}

  /** @name Numerical Methods (Friends)
   */
  //@{

  /// returns the sum of the diagonal elements A(i,i) over all i
  friend Real trace(const Sparsemat& A);               //=sum(diag(A))
  ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  friend Real ip(const Sparsemat& A, const Sparsemat& B); //=trace(B^t*A)
  ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  friend Real ip(const Sparsemat& A, const Matrix& B); //=trace(B^t*A)
  ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  friend inline Real ip(const Matrix& A, const Sparsemat& B);

  ///returns a row vector holding the sum over all rows, i.e., (1 1 ... 1)*A
  friend Matrix sumrows(const Sparsemat& A);   //=(1 1 1 ... 1)*A
  ///returns a column vector holding the sum over all columns, i.e., A*(1 1 ... 1)^T
  friend Matrix sumcols(const Sparsemat& A);   //=A*(1 1 ... 1)^t
  ///returns the sum over all elements of A, i.e., (1 1 ... 1)*A*(1 1 ... 1)^T
  friend inline Real sum(const Sparsemat& A);  //=(1 1 ... 1)*A*(1 1 ... 1)^t

  //@}    

  //---------------------------------------------
  //----  Comparisons / Max / Min / sort / find
  //---------------------------------------------

  /** @name Equal (Members)
   */
  //@{

  /// returns 1 if both matrices are identical, 0 otherwise
  friend int equal(const Sparsemat& A, const Sparsemat& B,Real eqtol);

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

  ///output format: nr nc nz \\n i1 j1 val1\\n i2 j2 val2\\n ... inz jnz valnz\\n 
    friend std::ostream& operator<<(std::ostream& o,const Sparsemat &v);
  ///input format: nr nc nz \\n i1 j1 val1\\n i2 j2 val2\\n ... inz jnz valnz\\n 
    friend std::istream& operator>>(std::istream& i,Sparsemat &v);

  //@}

};

//@}


// **************************************************************************
//                make non inline friends available outside
// **************************************************************************

  void swap(Sparsemat& A, Sparsemat& B); 
  Sparsemat& xbpeya(Sparsemat& x,const Sparsemat& y,Real alpha,Real beta,int ytrans);
  Matrix& genmult(const Sparsemat& A,const Matrix& B,Matrix &C,
			 Real alpha,Real beta, int atrans,int btrans);
  Matrix& genmult(const Matrix& A,const Sparsemat& B,Matrix &C,
			 Real alpha,Real beta, int atrans,int btrans);
  Matrix& genmult(const Sparsemat& A,const Sparsemat& B,Matrix &C,
			 Real alpha,Real beta, int atrans,int btrans);
  Sparsemat operator*(const Sparsemat& A,const Sparsemat& B);
  Matrix& genmult(const Symmatrix& A,const Sparsemat& B,Matrix& C,
			 Real alpha,Real beta,int btrans);
  Matrix& genmult(const Sparsemat& A,const Symmatrix& B,Matrix& C,
			 Real alpha,Real beta, int atrans);
  Matrix& genmult(const Sparsesym& A,const Sparsemat& B,Matrix &C,
			 Real alpha,Real beta, int btrans);
  Matrix& genmult(const Sparsemat& A,const Sparsesym& B,Matrix &C,
			 Real alpha,Real beta,int atrans);
  Symmatrix& rankadd(const Sparsemat& A,Symmatrix& C,
			    Real alpha,Real beta,int trans);
  Symmatrix& rank2add(const Sparsemat& A,const Matrix& B,Symmatrix& C,
			     Real alpha,Real beta,int trans);
  Sparsemat abs(const Sparsemat& A);
  Real trace(const Sparsemat& A);               //=sum(diag(A))
  Real ip(const Sparsemat& A, const Sparsemat& B); //=trace(B^t*A)
  Real ip(const Sparsemat& A, const Matrix& B); //=trace(B^t*A)
  Matrix sumrows(const Sparsemat& A);   //=(1 1 1 ... 1)*A
  Matrix sumcols(const Sparsemat& A);   //=A*(1 1 ... 1)^t  
  std::ostream& operator<<(std::ostream& o,const Sparsemat &v);
  std::istream& operator>>(std::istream& i,Sparsemat &v);


// **************************************************************************
//                   implementation of inline functions
// **************************************************************************

inline void Sparsemat::init_to_zero()
{
 mtype=MTsparse;
 nr=nc=0;
 tol=SPARSE_ZERO_TOL;
#if (CONICBUNDLE_DEBUG>=1)
 is_init=true;
#endif
}

//initialize

inline Sparsemat& Sparsemat::init(const Sparsemat& A,Real d)
{ return xeya(A,d);}
inline Sparsemat& Sparsemat::init(const Matrix& A,Real d)
{ return xeya(A,d);}
inline Sparsemat& Sparsemat::init(const Indexmatrix& A,Real d)
{ return xeya(A,d);}


inline Sparsemat& Sparsemat::init(Integer r,Integer c)
{
chk_range(r,c,-1,-1); 
nr=r; nc=c;
colinfo.init(0,0,Integer(0)); colindex.init(0,0,Integer(0)); colval.init(0,0,0.);
rowinfo.init(0,0,Integer(0)); rowindex.init(0,0,Integer(0)); rowval.init(0,0,0.);
chk_set_init(*this,1);
return *this;
}

inline Sparsemat::Sparsemat()
{ init_to_zero(); chk_set_init(*this,1);}

  inline Sparsemat::Sparsemat(const Sparsemat& A,double d):Memarrayuser()
{ init_to_zero(); xeya(A,d);}
inline Sparsemat::Sparsemat(const Matrix& A,double d)
{ init_to_zero(); xeya(A,d);}
inline Sparsemat::Sparsemat(const Indexmatrix& A,double d)
{ init_to_zero(); xeya(A,d);}

inline Sparsemat::Sparsemat(Integer r,Integer c)
{ init_to_zero(); init(r,c); }

inline Sparsemat::Sparsemat(Integer in_nr,Integer in_nc,Integer nz,
       const Integer *ini,const Integer *inj,const Real* va)
{ init_to_zero(); init(in_nr,in_nc,nz,ini,inj,va);}

inline Sparsemat::Sparsemat(Integer in_nr,Integer in_nc,Integer nz,
              const Indexmatrix& ini,const Indexmatrix& inj, const Matrix& va)
{ init_to_zero(); init(in_nr,in_nc,nz,ini,inj,va);}

inline Real Sparsemat::operator()(Integer i) const
{return (*this)(i%nr,i/nr);}           
inline Real Sparsemat::operator[](Integer i) const
{return (*this)(i%nr,i/nr);}
    
inline Sparsemat& Sparsemat::transpose()
{chk_init(*this);Integer h=nr; nr=nc;nc=h;swap(colinfo,rowinfo);
 swap(colindex,rowindex);swap(colval,rowval);return *this;}          
    
inline Sparsemat concat_right(const Sparsemat& A,const Sparsemat& B)
{Sparsemat C(A);C.concat_right(B);return C;}
inline Sparsemat concat_below(const Sparsemat& A,const Sparsemat& B)
{Sparsemat C(A);C.concat_below(B);return C;}


inline Sparsemat& Sparsemat::operator=(const Sparsemat &A)
{ return xeya(A); }
inline Sparsemat& Sparsemat::operator+=(const Sparsemat &A)
{ return xbpeya(*this,A,1.,1.); }
inline Sparsemat& Sparsemat::operator-=(const Sparsemat &A)
{ return xbpeya(*this,A,-1.,1.); }
inline Sparsemat& Sparsemat::operator*=(Real d)
{chk_init(*this); colval*=d; rowval*=d; return *this;}
inline Sparsemat& Sparsemat::operator/=(Real d)
{chk_init(*this);colval/=d;rowval/=d;return *this;}
inline Sparsemat Sparsemat::operator-() const
{return Sparsemat(*this,-1.); }
inline Sparsemat operator+(const Sparsemat &A,const Sparsemat& B)
{ Sparsemat C(A); return xbpeya(C,B,1.,1.); }
inline Sparsemat operator-(const Sparsemat &A,const Sparsemat& B)
{ Sparsemat C(A); return xbpeya(C,B,-1.,1.); }
inline Sparsemat operator*(const Sparsemat& A,Real d)
{return Sparsemat(A,d);}
inline Sparsemat operator*(Real d,const Sparsemat& A) 
{return Sparsemat(A,d);}
inline Sparsemat operator/(const Sparsemat& A,Real d) 
{return Sparsemat(A,1./d);}

inline Sparsemat& Sparsemat::operator=(const Matrix &A)
{return xeya(A);}
inline Sparsemat& Sparsemat::operator=(const Indexmatrix &A)
{return xeya(A);}


inline Matrix operator*(const Sparsemat& A,const Matrix& B) 
    {Matrix C;return genmult(A,B,C);}
inline Matrix operator*(const Matrix& A,const Sparsemat& B) 
    {Matrix C;return genmult(A,B,C);}

inline Matrix operator+(const Sparsemat& A,const Matrix& B)
    {Matrix C(B);return C.xpeya(A);}
inline Matrix operator+(const Matrix& A,const Sparsemat& B)
    {Matrix C(A);return C.xpeya(B);}
inline Matrix operator-(const Sparsemat& A,const Matrix& B)
    {Matrix C(B,-1.);return C.xpeya(A);}
inline Matrix operator-(const Matrix& A,const Sparsemat& B)
    {Matrix C(A);return C.xpeya(B,-1.);}


inline Sparsemat::Sparsemat(const Symmatrix& A,Real d)
{ init_to_zero(); xeya(Matrix(A),d);}
inline Sparsemat& Sparsemat::init(const Symmatrix &A,Real d)
{return xeya(Matrix(A),d);}
inline Sparsemat& Sparsemat::operator=(const Symmatrix &A)
{return xeya(Matrix(A));}

inline Real ip(const Matrix& A, const Sparsemat& B)
{return ip(B,A);}

inline Real sum(const Sparsemat& A)
{return sum(A.colval);}

inline Sparsemat transpose(const Sparsemat& A)
{Sparsemat B(A); B.transpose(); return B;}



inline Matrix::Matrix(const Sparsemat& A, Real d)        
{init_to_zero(); xeya(A,d);}

inline Matrix& Matrix::init(const Sparsemat& A, Real d)  
{return xeya(A,d);}

inline Matrix& Matrix::operator=(const Sparsemat & A)       
{return xeya(A);}

inline Matrix& Matrix::operator*=(const Sparsemat &A)
{ Matrix C; return *this=genmult(*this,A,C);}

inline Matrix& Matrix::operator+=(const Sparsemat &A)       
{return xpeya(A);}

inline Matrix& Matrix::operator-=(const Sparsemat &A)       
{return xpeya(A,-1.);}


  /** @name Equal (Members)
   */
  //@{

  /// returns 1 if both matrices are identical, 0 otherwise
  int equal(const Sparsemat& A, const Sparsemat& B,Real eqtol=1e-10);

  //@}    


}

//in order to define the inlines associated with the other matrix
//classes, the hxx-files of these are included here, as well. 

#ifndef CH_MATRIX_CLASSES__SPARSSYM_HXX
#include "sparssym.hxx"
#endif

#endif

