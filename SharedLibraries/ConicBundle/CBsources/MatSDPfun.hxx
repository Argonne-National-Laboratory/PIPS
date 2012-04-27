/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatSDPfun.hxx

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



#ifndef CONICBUNDLE_MATSDPFUN_HXX
#define CONICBUNDLE_MATSDPFUN_HXX

//------------------------------------------------------------

#include <map>
#include "BaseSDPOracle.hxx"
#include "matrix.hxx"
#include "sparssym.hxx"
#include "coeffmat.hxx"
#include "bigmat.hxx"

//------------------------------------------------------------

/**@brief   oracle for SDP with matrix classes
	@author  C. Helmberg
*/	
namespace ConicBundle {

  class DenseSDPPrimal: public SDPPrimal, public CH_Matrix_Classes::Symmatrix
  {
  private:

  public:
    DenseSDPPrimal(){}
    DenseSDPPrimal(const DenseSDPPrimal& symmat) : SDPPrimal(), CH_Matrix_Classes::Symmatrix(symmat) {}
    DenseSDPPrimal(const CH_Matrix_Classes::Symmatrix& symmat) : SDPPrimal(), CH_Matrix_Classes::Symmatrix(symmat) {}
    const DenseSDPPrimal& operator=(const CH_Matrix_Classes::Symmatrix& symmat)
    { init(symmat); return *this; }
 
     /// returns a newly generated identical Object
    virtual PrimalData* clone_primal_data()
    { return new DenseSDPPrimal(*this); }

    /// this assign is only feasible for selected derivations of PrimalData 
    virtual int assign_primal_data(const PrimalData& pd)
    {
      const DenseSDPPrimal* p=dynamic_cast<const DenseSDPPrimal*>(&pd);
      assert(p!=0);
      *this=*p;
      return 0;
    }
      

    /// set aggregate=P*P^T
    virtual int assign_Gram_matrix(const CH_Matrix_Classes::Matrix& P)
    { CH_Matrix_Classes::rankadd(P,*this); return 0; }
    
    /** multiply this with myfactor and add itsfactor*it to this
        (it must also be a DenseSDPPrimal) */
    virtual int aggregate_primal_data(double myfactor,double itsfactor,const PrimalData& it)
    { 
      const DenseSDPPrimal* pd=dynamic_cast<const DenseSDPPrimal*>(&it);
      assert(pd!=0);
      CH_Matrix_Classes::xbpeya(*this,*pd,itsfactor,myfactor);
      return 0;
    }

    /// multiply this with myfactor and add itsfactor*P*P^T to this
    virtual int aggregate_Gram_matrix(double myfactor,double itsfactor,const CH_Matrix_Classes::Matrix& P)
    { CH_Matrix_Classes::rankadd(P,*this,itsfactor,myfactor); return 0; }

  };


  class SparseSDPPrimal: public SDPPrimal, public CH_Matrix_Classes::Sparsesym
  {
  public:
    SparseSDPPrimal(const CH_Matrix_Classes::Sparsesym& sps) : SDPPrimal(), CH_Matrix_Classes::Sparsesym(sps) {}
    SparseSDPPrimal(const SparseSDPPrimal& pr) : SDPPrimal(), CH_Matrix_Classes::Sparsesym(pr)  {} 
    ~SparseSDPPrimal(){}
    SparseSDPPrimal& operator=(const CH_Matrix_Classes::Sparsesym& sdp)  
    { init(sdp); return *this; }
    
 
 
     /// returns a newly generated identical Object
    PrimalData* clone_primal_data()
    { return new SparseSDPPrimal(*this); }

    /// this assign is only feasible for selected derivations of PrimalData 
    int assign_primal_data(const PrimalData& pd)
    { 
      const SparseSDPPrimal* p=dynamic_cast<const SparseSDPPrimal*>(&pd);
      assert(p!=0);
      *this = *p;
      return 0;
    }

    /// for each element aij in the support set aij=<P.row(i),P.row(j)> 
    int assign_Gram_matrix(const CH_Matrix_Classes::Matrix& P)
    { CH_Matrix_Classes::support_rankadd(P,*this); return 0; }
    
    /// multiply this with myfactor and add itsfactor*it to this
    int aggregate_primal_data(double myfactor,double itsfactor,const PrimalData& it)
    {
      const SparseSDPPrimal* pd=dynamic_cast<const SparseSDPPrimal*>(&it);
      assert(pd!=0);
      CH_Matrix_Classes::xbpeya(*this,*pd,itsfactor,myfactor);
      return 0;
    }

    /// multiply this with myfactor and add itsfactor*P*P^T to this
    int aggregate_Gram_matrix(double myfactor,double itsfactor,const CH_Matrix_Classes::Matrix& P)
    { CH_Matrix_Classes::support_rankadd(P,*this,itsfactor,myfactor); return 0; }

  };
  

  class GramSparseSDPPrimal: public SDPPrimal, public CH_Matrix_Classes::Sparsesym
  {
   protected:
    CH_Matrix_Classes::Matrix grammatrix;

   public:
    GramSparseSDPPrimal(const CH_Matrix_Classes::Sparsesym& sps) : SDPPrimal(), CH_Matrix_Classes::Sparsesym(sps) {}
    GramSparseSDPPrimal(const GramSparseSDPPrimal& pr) : 
      SDPPrimal(), 
      CH_Matrix_Classes::Sparsesym(pr),
      grammatrix(pr.grammatrix)  {} 
    ~GramSparseSDPPrimal(){}
    GramSparseSDPPrimal& operator=(const CH_Matrix_Classes::Sparsesym& sdp)  
    { init(sdp); grammatrix.init(0,0,0.); return *this; }
    GramSparseSDPPrimal& operator=(const GramSparseSDPPrimal& sdp)  
    { init(sdp); grammatrix=sdp.grammatrix; return *this; }
    
 
    const CH_Matrix_Classes::Matrix& get_grammatrix() const { return  grammatrix; }
 
     /// returns a newly generated identical Object
    PrimalData* clone_primal_data()
    { return new GramSparseSDPPrimal(*this); }

    /// this assign is only feasible for selected derivations of PrimalData 
    int assign_primal_data(const PrimalData& pd)
    { 
      const GramSparseSDPPrimal* p=dynamic_cast<const GramSparseSDPPrimal*>(&pd);
      if(p!=0){
	*this = *p;
	return 0;
      }
      const SparseSDPPrimal* sp=dynamic_cast<const SparseSDPPrimal*>(&pd);
      if(sp!=0){
	*this = dynamic_cast<const CH_Matrix_Classes::Sparsesym&>(*sp);
	return 0;
      }
      return 1;
    }

    /// for each element aij in the support set aij=<P.row(i),P.row(j)> 
    int assign_Gram_matrix(const CH_Matrix_Classes::Matrix& P)
    { grammatrix=P; CH_Matrix_Classes::support_rankadd(CH_Matrix_Classes::Matrix(grammatrix.rowdim(),1,0.),*this,0.); return 0; }
    
    /// multiply this with myfactor and add itsfactor*it to this
    int aggregate_primal_data(double myfactor,double itsfactor,const PrimalData& it)
    {
      const GramSparseSDPPrimal* pd=dynamic_cast<const GramSparseSDPPrimal*>(&it);
      if (pd!=0){
	if ((grammatrix.dim()!=0)&&(myfactor!=1.)) {
	  grammatrix*=sqrt(myfactor);
	}
	CH_Matrix_Classes::xbpeya(*this,*pd,itsfactor,myfactor);
	if (pd->grammatrix.dim()!=0) CH_Matrix_Classes::support_rankadd(pd->grammatrix,*this,itsfactor);
	return 0;
      }
      const SparseSDPPrimal* ps=dynamic_cast<const SparseSDPPrimal*>(&it);
      if (ps!=0){
	if ((grammatrix.dim()!=0)&&(myfactor!=1.)) {
	  grammatrix*=sqrt(myfactor);
	}
	CH_Matrix_Classes::xbpeya(*this,*ps,itsfactor,myfactor);
	return 0;
      }
      return 1;
    }

    /// multiply this with myfactor and add itsfactor*P*P^T to this
    int aggregate_Gram_matrix(double myfactor,double itsfactor,const CH_Matrix_Classes::Matrix& P)
    { 
      if (grammatrix.dim()!=0) {
	CH_Matrix_Classes::support_rankadd(grammatrix,*this);
	grammatrix.init(0,0,0.);
      }
      CH_Matrix_Classes::support_rankadd(P,*this,itsfactor,myfactor); 
      return 0; 
    }

  };
  


  class BlockSDPPrimal: public SDPPrimal
  {
  private:
    CH_Matrix_Classes::Indexmatrix Xdim;
    std::map<CH_Matrix_Classes::Integer,SDPPrimal*> primal;
  public:
    BlockSDPPrimal(const BlockSDPPrimal& pr);
    /** in this version control over the objects pointed to in pr 
	is passed to this. this will delete them upon its destruction */
    BlockSDPPrimal(const CH_Matrix_Classes::Indexmatrix& Xd,const std::map<CH_Matrix_Classes::Integer,SDPPrimal*>& pr);
 
    ~BlockSDPPrimal();
 
     /// returns a newly generated identical Object
    virtual PrimalData* clone_primal_data();

    /// this assign is only feasible for selected derivations of PrimalData 
    virtual int assign_primal_data(const PrimalData& pd);

    /// for each element aij in the support set aij=<P.row(i),P.row(j)> 
    virtual int assign_Gram_matrix(const CH_Matrix_Classes::Matrix& P);
    
    /** multiply this with myfactor and add itsfactor*it to this
        (it must also be a SparseSDPPrimal and on the same support) */
    virtual int aggregate_primal_data(double myfactor,double itsfactor,const PrimalData& it);

    /// multiply this with myfactor and add itsfactor*P*P^T to this
    virtual int aggregate_Gram_matrix(double myfactor,double itsfactor,const CH_Matrix_Classes::Matrix& P);

    /// get dimension information
    CH_Matrix_Classes::Integer get_nblocks() const {return Xdim.dim();}
    /// get dimension information
    CH_Matrix_Classes::Integer block_dim(CH_Matrix_Classes::Integer i) const {return Xdim(i);}
    /// get aggregate information
    SDPPrimal* block(CH_Matrix_Classes::Integer i) const;
  };


  class MaxEigSolver;

  typedef std::map<CH_Matrix_Classes::Integer,Coeffmat*> SparseSDPCoeffVector;
  typedef std::map<CH_Matrix_Classes::Integer,SparseSDPCoeffVector*> SparseSDPCoeffMatrix;

  class MatrixSDPfunction: public BaseSDPOracle 
  {
  protected:
    CH_Matrix_Classes::Real trace_value;
    SDPtrace constant_trace;
    SparseSDPCoeffVector C;
    CH_Matrix_Classes::Indexmatrix Xdim;                         //dimension of each variable

    //one SparseSDPCoeffVector for each nonzero row
    SparseSDPCoeffMatrix opA_rowrep;  

    //one SparseSDPCoeffVector for each nonzero column
    SparseSDPCoeffMatrix opA_colrep;  

    CH_Matrix_Classes::Matrix b;

    //for each variable: number of dense coefficient matrices in C and opA
    CH_Matrix_Classes::Indexmatrix n_dense_coeffmat;  

    /** This SDPPrimal can be set from outside and serves for
        generating further primals by cloning */
    SDPPrimal* generating_primal; 
    
    std::vector<MaxEigSolver*> maxeigsolver;
    CH_Matrix_Classes::Matrix last_bigmat_y; 

    CH_Matrix_Classes::Integer maxvecs;

    CH_Matrix_Classes::Matrix tmpmat;

    std::ostream* out;
    CH_Matrix_Classes::Integer print_level;
    
  public:

    MatrixSDPfunction();
    ~MatrixSDPfunction();

    void set_maxvecs(CH_Matrix_Classes::Integer i){maxvecs=CH_Matrix_Classes::max(i,1);}


    //----------- Oracle Implementation of BaseSDPOracle ----------

    void get_trace_constraint(CH_Matrix_Classes::Real& trace_val,SDPtrace& trace_stat);
 
    void adjust_multiplier(CH_Matrix_Classes::Real& mult)
    { 
      if ((constant_trace==SDPtrace_unbounded)&&(mult>0.))
	trace_value=mult;
      else mult=trace_value;
    }

    int gramip(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Real& ip_C,CH_Matrix_Classes::Matrix& ip_opA);

    int gramip_opA(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Matrix& vec);

    const CH_Matrix_Classes::Matrix& rhs() const;
        
    int project_C(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Symmatrix& S);

    int project(const CH_Matrix_Classes::Integer i, const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Symmatrix& S);

    SDPPrimal* init_primal(const CH_Matrix_Classes::Matrix& P);
    
    int primalip_opA(const SDPPrimal* primal,const CH_Matrix_Classes::Indexmatrix& ind,
			     CH_Matrix_Classes::Matrix& ip_opA);
  
    //in addtion, because it is useful
    int primalip_C(const SDPPrimal* primal,CH_Matrix_Classes::Real& value);
  
    
    int evaluate(const CH_Matrix_Classes::Matrix& current_point, const CH_Matrix_Classes::Matrix& bundlevecs, 
		 const double relprec,const double Ritz_bound,
		 CH_Matrix_Classes::Matrix& Ritz_vectors,CH_Matrix_Classes::Matrix&  Ritz_values);

    int evaluate_projection(const CH_Matrix_Classes::Matrix& current_point, const CH_Matrix_Classes::Matrix& P, 
			     const double relprec, const double Ritz_bound,
			     CH_Matrix_Classes::Matrix& projected_Ritz_vectors,
			     CH_Matrix_Classes::Matrix& projected_Ritz_values);

    int compute_projection(const CH_Matrix_Classes::Matrix& current_point, const CH_Matrix_Classes::Matrix& P,  CH_Matrix_Classes::Symmatrix &S);

    int eig(const CH_Matrix_Classes::Matrix& current_point,CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Matrix& d);

    //------------------  routines for external computations ----------------
   
    /// computes all eigenvectors P and eigenvalues d of block ind of C-opAt(y) 
    int eig(const CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Integer ind,
	    CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Matrix& d) const; 

    /// returns a pointer to the internal Bigmatrix representation of block ind of C-opAt(y) 
    const Bigmatrix& get_bigmat(const CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Integer ind);

    //------------------  routines for modifying the problem ----------------

    /// set the basic constraint <I,X> <=/= trace; #trace# must be nonnegative
    int set_trace(const CH_Matrix_Classes::Real trace, SDPtrace constant_trace);

   
    /** defines in what form primal matrices should be aggregated.
        If the argument is NULL then no primal aggregation will take place.
	The control over the generating primal is
        passed over to this. This will delete an existing generating primal 
        whenever a new generating primal is set or upon destruction. */
    int set_generating_primal(SDPPrimal* gen_primal);

    const SDPPrimal* get_generating_primal(void)
    {return generating_primal;}
        
    /** add new variables and their coefficients for existing constraints.
        If #append_Xdim# has dimension k then k matrix variables are added. 
	#append_Xdim(k)# gives the order of matrix k, #append_C.find(k)->second# 
	points to the cost matrix and #opA_columns.find(k)->second# 
	gives the map for the constraint cofficient matrices of matrix 
	variable k. It is not possible to introduce new constraints by 
	this routine. 
        The control over the coefficient matrices pointed to by C and opA
        is passed over to #*this#. #*this# will delete the coefficient matrices
	if replaced by a new one or on destruction of #*this#. */
    int append_variables(const CH_Matrix_Classes::Indexmatrix append_Xdim,
			 const SparseSDPCoeffVector& append_C,
			 const std::map<CH_Matrix_Classes::Integer,SparseSDPCoeffVector>& opA_columns);


    /** append constraints to the problem. The #SparseSDPCoeffVector#
	of #append_opA[i]# belongs to matrix variable i. 
        The control over all coefficient matrices pointed to in 
	#append_opA# is passed over to this. This will delete them
        when delete constraints is called or on destruction.
        The other values are copied. */
    int append_constraints(const std::map<CH_Matrix_Classes::Integer,SparseSDPCoeffVector>& opA_rows,
			   const CH_Matrix_Classes::Matrix& append_b);

    /** The new assignment of the remaining 
	variables is described by #map_to_old# so that aftwards
	position i is assigned to the variable previously indexed
        by #map_to_old(i)#. Generating copies is allowed.  
        If #deleted_cols!=0# or #deleted_C!=0# or #deleted_cols!=0#
	appropriate deleted column data is stored as new
        columns in the respective variable pointed to (these are
        assumed to be empty on input), otherwise it is deleted. */        
    int reassign_variables(const CH_Matrix_Classes::Indexmatrix& map_to_old,
			   CH_Matrix_Classes::Indexmatrix* deleted_Xdim=0,
			   SparseSDPCoeffVector* deleted_C=0,
			   SparseSDPCoeffMatrix* deleted_cols=0);
 
    /** delete variables. If #map_to_old# is not null then
        the correspondence of the new indices to the old indices of
        the remaining variables is described in #map_to_old# so 
        that now position i corresponds to the variable previously 
        indexed by #map_to_old(i)#.   
        If #deleted_cols!=0# or #deleted_C!=0# or #deleted_cols!=0#
	appropriate deleted column data is stored as new
        columns in the respective variable pointed to (these are
        assumed to be empty on input), otherwise it is deleted. */   
    int delete_variables(const CH_Matrix_Classes::Indexmatrix& delete_indices,
			 CH_Matrix_Classes::Indexmatrix* map_to_old=0,
			 CH_Matrix_Classes::Indexmatrix* deleted_Xdim=0,
			 SparseSDPCoeffVector* deleted_C=0,
			 SparseSDPCoeffMatrix* deleted_cols=0);
 
    /** The new assignment of the remaining 
	constraints is described by #map_to_old# so that aftwards
	position i is assigned the constraint previously indexed
        by #map_to_old(i)#. Generating copies is allowed. 
        If #deleted_rows!=0# or #deleted_rhs!=0# the appropriate
        deleted row data is appended as new rows in the respective 
	variable pointed to (these are assumed to be empty on input), 
	otherwise it is deleted. */        
    int reassign_constraints(const CH_Matrix_Classes::Indexmatrix& map_to_old,
			     SparseSDPCoeffMatrix* deleted_rows=0,
			     CH_Matrix_Classes::Matrix* deleted_rhs=0);
 
    /** delete constraints. If #map_to_old# is not null then
        the correspondence of the new indices to the old indices of
        the remaining constraints is described in #map_to_old# so 
        that now position i corresponds to the constraint previously 
        indexed by #map_to_old(i)#. 
        If #deleted_rows!=0# or #deleted_rhs!=0# the appropriate
        deleted row data is appended as new rows in the respective 
	variable pointed to (these are assumed to be empty on input), 
	otherwise it is deleted. */        
    int delete_constraints(const CH_Matrix_Classes::Indexmatrix& delete_indices,
			   CH_Matrix_Classes::Indexmatrix* map_to_old=0,
			   SparseSDPCoeffMatrix* deleted_rows=0,
			   CH_Matrix_Classes::Matrix* deleted_rhs=0);
 
   
    //------------------  routines for querying the problem ----------------

    const CH_Matrix_Classes::Indexmatrix& get_Xdim() const {return Xdim;}
    
    /** returns the pointer to the coefficient matrix corresponding
     to constraint constr_nr and matrix variable block block_nr.
     The returned pointer may be NULL if the coefficient matrix
     is not specified (i.e. if it is regarded as the zero matrix)
     or the indices are out of range.
     The coefficient matrices of the objective are obtained for
     constr_nr==-1.
    */
    const Coeffmat* get_coeffmat(CH_Matrix_Classes::Integer constr_nr,CH_Matrix_Classes::Integer block_nr) const;

    /** returns the row representation of the coefficient matrices
     (each entry of the map represents a row by a SparseSDPCoeffVector).
    */
    const SparseSDPCoeffMatrix& get_opA_rowrep() const {return opA_rowrep;}
    
    /** returns the column representation of the coefficient matrices
     (each entry of the map represents a column by a SparseSDPCoeffVector).
    */
    const SparseSDPCoeffMatrix& get_opA_colrep() const {return opA_colrep;}
    

    //------------------  routines for IO ----------------

    void  set_out(std::ostream* o=0,int pril=1);

    std::ostream& print_problem_data(std::ostream& out);
    std::istream& read_problem_data(std::istream& in);

  };
  //@}

  int convert_row(const CH_Matrix_Classes::Indexmatrix& sdpdim,const CH_Matrix_Classes::Sparsemat& row,SparseSDPCoeffVector&);
  int convert_to_rows(const CH_Matrix_Classes::Indexmatrix& sdpdim,const CH_Matrix_Classes::Sparsemat& A,std::map<CH_Matrix_Classes::Integer,SparseSDPCoeffVector>& opA_rows);
  int convert_to_cols(const CH_Matrix_Classes::Indexmatrix& sdpdim,const CH_Matrix_Classes::Sparsemat& A,std::map<CH_Matrix_Classes::Integer,SparseSDPCoeffVector>& opA_cols);
}

#endif

