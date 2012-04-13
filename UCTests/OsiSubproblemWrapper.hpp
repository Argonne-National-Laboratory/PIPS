#ifndef OSISUBPROBLEMWRAPPER_HPP
#define OSISUBPROBLEMWRAPPER_HPP

#include "OsiSolverInterface.hpp"
#include "stochasticInput.hpp"
#include <vector>
#include "CoinFinite.hpp"

// wrap a subproblem with 1st stage and 1 scenario so we can give it to Cgl to generate cuts
class OsiSubproblemWrapper : public OsiSolverInterface {

public:

	void setCurrentSolution(const std::vector<double> stg1, const std::vector<double> stg2);
	void setCurrentReducedCosts(const std::vector<double> stg1, const std::vector<double> stg2);

  //---------------------------------------------------------------------------

  ///@name Solve methods 
  //@{
    /// Solve initial LP relaxation 
    virtual void initialSolve() {assert(0);}

    /*! \brief Resolve an LP relaxation after problem modification

      Note the `re-' in `resolve'. initialSolve() should be used to solve the
      problem for the first time.
    */
    virtual void resolve() {assert(0);}

    /// Invoke solver's built-in enumeration algorithm
    virtual void branchAndBound() { assert(0);}

   
    //@}

  //---------------------------------------------------------------------------
  ///@name Methods returning info on how the solution process terminated
  //@{
    /// Are there numerical difficulties?
    virtual bool isAbandoned() const {assert(0); return false;}
    /// Is optimality proven?
    virtual bool isProvenOptimal() const {assert(0); return false;}
    /// Is primal infeasibility proven?
    virtual bool isProvenPrimalInfeasible() const {assert(0); return false;}
    /// Is dual infeasibility proven?
    virtual bool isProvenDualInfeasible() const { assert(0); return false;}
    /// Iteration limit reached?
    virtual bool isIterationLimitReached() const { assert(0); return false;}
  //@}

  //---------------------------------------------------------------------------
  /** \name Warm start methods

    Note that the warm start methods return a generic CoinWarmStart object.
    The precise characteristics of this object are solver-dependent. Clients
    who wish to maintain a maximum degree of solver independence should take
    care to avoid unnecessary assumptions about the properties of a warm start
    object.
  */
  //@{
    /*! \brief Get an empty warm start object
      
      This routine returns an empty warm start object. Its purpose is
      to provide a way for a client to acquire a warm start object of the
      appropriate type for the solver, which can then be resized and modified
      as desired.
    */

    virtual CoinWarmStart *getEmptyWarmStart () const {assert(0); return 0;}

    /** \brief Get warm start information.

      Return warm start information for the current state of the solver
      interface. If there is no valid warm start information, an empty warm
      start object wil be returned.
    */
    virtual CoinWarmStart* getWarmStart() const {assert(0); return 0;}
    /** \brief Set warm start information.
    
      Return true or false depending on whether the warm start information was
      accepted or not.
      By definition, a call to setWarmStart with a null parameter should
      cause the solver interface to refresh its warm start information
      from the underlying solver.
   */
    virtual bool setWarmStart(const CoinWarmStart* warmstart) {assert(0); return false;}
  //@}


  //---------------------------------------------------------------------------
  /**@name Problem query methods

   Querying a problem that has no data associated with it will result in
   zeros for the number of rows and columns, and NULL pointers from the
   methods that return vectors.

   Const pointers returned from any data-query method are valid as long as
   the data is unchanged and the solver is not called.
  */
  //@{
    /// Get the number of columns
    virtual int getNumCols() const { return nvar1 + nvar2;}

    /// Get the number of rows
    virtual int getNumRows() const { return ncons2; }

    /// Get the number of nonzero elements
    virtual int getNumElements() const {assert(0); return 0;}


    /// Get a pointer to an array[getNumCols()] of column lower bounds
    virtual const double * getColLower() const { return &collb[0];} 

    /// Get a pointer to an array[getNumCols()] of column upper bounds
    virtual const double * getColUpper() const { return &colub[0];} 

    /*! \brief Get a pointer to an array[getNumRows()] of row constraint senses.

      <ul>
      <li>'L': <= constraint
      <li>'E': =  constraint
      <li>'G': >= constraint
      <li>'R': ranged constraint
      <li>'N': free constraint
      </ul>
    */
    virtual const char * getRowSense() const { return &rowSense[0];} 

    /*! \brief Get a pointer to an array[getNumRows()] of row right-hand sides

      <ul>
	<li> if getRowSense()[i] == 'L' then
	     getRightHandSide()[i] == getRowUpper()[i]
	<li> if getRowSense()[i] == 'G' then
	     getRightHandSide()[i] == getRowLower()[i]
	<li> if getRowSense()[i] == 'R' then
	     getRightHandSide()[i] == getRowUpper()[i]
	<li> if getRowSense()[i] == 'N' then
	     getRightHandSide()[i] == 0.0
      </ul>
    */
    virtual const double * getRightHandSide() const { return &rhs[0]; } // TODO

    /*! \brief Get a pointer to an array[getNumRows()] of row ranges.

      <ul>
	  <li> if getRowSense()[i] == 'R' then
		  getRowRange()[i] == getRowUpper()[i] - getRowLower()[i]
	  <li> if getRowSense()[i] != 'R' then
		  getRowRange()[i] is 0.0
	</ul>
    */
    virtual const double * getRowRange() const {assert(0); return 0;} // TODO

    /// Get a pointer to an array[getNumRows()] of row lower bounds
    virtual const double * getRowLower() const { return &rowlb[0];} 
    /// Get a pointer to an array[getNumRows()] of row upper bounds
    virtual const double * getRowUpper() const { return &rowub[0];} 

    /*! \brief Get a pointer to an array[getNumCols()] of objective
	       function coefficients.
    */
    virtual const double * getObjCoefficients() const {/*printf("warning: asked for objective coefficients\n");*/ return 0;}

    /*! \brief Get the objective function sense
    
      -  1 for minimisation (default)
      - -1 for maximisation
    */
    virtual double getObjSense() const { return 1;}

    /// Return true if the variable is continuous
    virtual bool isContinuous(int colIndex) const { return !isInteger[colIndex];} 



  
    /// Get a pointer to a row-wise copy of the matrix
    virtual const CoinPackedMatrix * getMatrixByRow() const { return &rMat;} 

    /// Get a pointer to a column-wise copy of the matrix
    virtual const CoinPackedMatrix * getMatrixByCol() const { return &cMat;} 


    /// Get the solver's value for infinity
    virtual double getInfinity() const { return COIN_DBL_MAX; }
  //@}
    
  /**@name Solution query methods */
  //@{
    /// Get a pointer to an array[getNumCols()] of primal variable values
    virtual const double * getColSolution() const {assert(colsol.size()); return &colsol[0];}

  

    /// Get pointer to array[getNumRows()] of dual variable values
    virtual const double * getRowPrice() const {assert(0); return 0; }

    /// Get a pointer to an array[getNumCols()] of reduced costs
    virtual const double * getReducedCost() const {assert(coldualsol.size()); return &coldualsol[0];}

    /** Get a pointer to array[getNumRows()] of row activity levels.
    
      The row activity for a row is the left-hand side evaluated at the
      current solution.
    */
    virtual const double * getRowActivity() const { assert(colsol.size()); return &rowActivity[0];} // TODO

    /// Get the objective function value.
    virtual double getObjValue() const {/*printf("warning, asked for objective value\n");*/ return 0;}

    /** Get the number of iterations it took to solve the problem (whatever
	`iteration' means to the solver).
    */
    virtual int getIterationCount() const {assert(0); return 0; }

    /** Get as many dual rays as the solver can provide. In case of proven
	primal infeasibility there should (with high probability) be at least
	one.

	The first getNumRows() ray components will always be associated with
	the row duals (as returned by getRowPrice()). If \c fullRay is true,
	the final getNumCols() entries will correspond to the ray components
	associated with the nonbasic variables. If the full ray is requested
	and the method cannot provide it, it will throw an exception.

	\note
	Implementors of solver interfaces note that the double pointers in
	the vector should point to arrays of length getNumRows() (fullRay =
	false) or (getNumRows()+getNumCols()) (fullRay = true) and they should
	be allocated with new[].

	\note
	Clients of solver interfaces note that it is the client's
	responsibility to free the double pointers in the vector using
	delete[]. Clients are reminded that a problem can be dual and primal
	infeasible.
    */
    virtual std::vector<double*> getDualRays(int maxNumRays,
					     bool fullRay = false) const {assert(0); return std::vector<double*>();}

    /** Get as many primal rays as the solver can provide. In case of proven
	dual infeasibility there should (with high probability) be at least
	one.
   
	\note
	Implementors of solver interfaces note that the double pointers in
	the vector should point to arrays of length getNumCols() and they
	should be allocated with new[].

	\note
	Clients of solver interfaces note that it is the client's
	responsibility to free the double pointers in the vector using
	delete[]. Clients are reminded that a problem can be dual and primal
	infeasible.
    */
    virtual std::vector<double*> getPrimalRays(int maxNumRays) const {assert(0); return std::vector<double*>();}

    
  //-------------------------------------------------------------------------
  /**@name Methods to modify the objective, bounds, and solution

     For functions which take a set of indices as parameters
     (\c setObjCoeffSet(), \c setColSetBounds(), \c setRowSetBounds(),
     \c setRowSetTypes()), the parameters follow the C++ STL iterator
     convention: \c indexFirst points to the first index in the
     set, and \c indexLast points to a position one past the last index
     in the set.
  
  */
  //@{
    /** Set an objective function coefficient */
    virtual void setObjCoeff( int elementIndex, double elementValue ) {assert(0);}

 

    /** Set the objective function sense.

        Use 1 for minimisation (default), -1 for maximisation.

	\note
	Implementors note that objective function sense is a parameter of
	the OSI, not a property of the problem. Objective sense can be
	set prior to problem load and should not be affected by loading a
	new problem.
    */
    virtual void setObjSense(double s) {assert(0);}
  

    /** Set a single column lower bound.
	Use -getInfinity() for -infinity. */
    virtual void setColLower( int elementIndex, double elementValue ) {assert(0);}
    
   
    /** Set a single column upper bound.
	Use getInfinity() for infinity. */
    virtual void setColUpper( int elementIndex, double elementValue ) {assert(0);}



    /** Set a single row lower bound.
	Use -getInfinity() for -infinity. */
    virtual void setRowLower( int elementIndex, double elementValue ) {assert(0);}
    
    /** Set a single row upper bound.
	Use getInfinity() for infinity. */
    virtual void setRowUpper( int elementIndex, double elementValue ) {assert(0);}


 
  
  
    /** Set the type of a single row */
    virtual void setRowType(int index, char sense, double rightHandSide,
			    double range) {assert(0);}


    /** Set the primal solution variable values

	colsol[getNumCols()] is an array of values for the primal variables.
	These values are copied to memory owned by the solver interface
	object or the solver.  They will be returned as the result of
	getColSolution() until changed by another call to setColSolution() or
	by a call to any solver routine.  Whether the solver makes use of the
	solution in any way is solver-dependent.
    */
    virtual void setColSolution(const double *colsol) {assert(0);}

    /** Set dual solution variable values

	rowprice[getNumRows()] is an array of values for the dual variables.
	These values are copied to memory owned by the solver interface
	object or the solver.  They will be returned as the result of
	getRowPrice() until changed by another call to setRowPrice() or by a
	call to any solver routine.  Whether the solver makes use of the
	solution in any way is solver-dependent.
    */
    virtual void setRowPrice(const double * rowprice) {assert(0);}

    //@}

  //-------------------------------------------------------------------------
  /**@name Methods to set variable type */
  //@{
    /** Set the index-th variable to be a continuous variable */
    virtual void setContinuous(int index) {assert(0);}
    /** Set the index-th variable to be an integer variable */
    virtual void setInteger(int index) {assert(0);}

  //@}
  //-------------------------------------------------------------------------
    
  //-------------------------------------------------------------------------
  /**@name Methods to modify the constraint system.

     Note that new columns are added as continuous variables.
  */
  //@{

    /** Add a column (primal variable) to the problem. */
    virtual void addCol(const CoinPackedVectorBase& vec,
			const double collb, const double colub,   
			const double obj) {assert(0);}

    /*! \brief Add a row (constraint) to the problem. */
    virtual void addRow(const CoinPackedVectorBase& vec,
			const double rowlb, const double rowub) {assert(0);}

   

    /*! \brief Add a row (constraint) to the problem. */
    virtual void addRow(const CoinPackedVectorBase& vec,
			const char rowsen, const double rowrhs,   
			const double rowrng) {assert(0);}
  


    /** \brief Remove a set of columns (primal variables) from the
	       problem.

      The solver interface for a basis-oriented solver will maintain valid
      warm start information if all deleted variables are nonbasic.
    */
    virtual void deleteCols(const int num, const int * colIndices) {assert(0);}
 


    /** \brief Delete a set of rows (constraints) from the problem.

      The solver interface for a basis-oriented solver will maintain valid
      warm start information if all deleted rows are loose.
    */
    virtual void deleteRows(const int num, const int * rowIndices) {assert(0);}

    virtual void writeMps (const char *filename,
			   const char *extension = "mps",
			   double objSense=0.0) const {assert(0);}

/** Apply a row cut (append to the constraint matrix). */
    virtual void applyRowCut( const OsiRowCut & rc ) {assert(0);}

    /** Apply a column cut (adjust the bounds of one or more variables). */
    virtual void applyColCut( const OsiColCut & cc ) {assert(0);}

//@}

  //---------------------------------------------------------------------------

  /**@name Methods for problem input and output */
  //@{
    /*! \brief Load in a problem by copying the arguments. The constraints on
	    the rows are given by lower and upper bounds.
	
	If a pointer is 0 then the following values are the default:
        <ul>
          <li> <code>colub</code>: all columns have upper bound infinity
          <li> <code>collb</code>: all columns have lower bound 0 
          <li> <code>rowub</code>: all rows have upper bound infinity
          <li> <code>rowlb</code>: all rows have lower bound -infinity
	  <li> <code>obj</code>: all variables have 0 objective coefficient
        </ul>

	Note that the default values for rowub and rowlb produce the
	constraint -infty <= ax <= infty. This is probably not what you want.
    */
    virtual void loadProblem (const CoinPackedMatrix& matrix,
			      const double* collb, const double* colub,   
			      const double* obj,
			      const double* rowlb, const double* rowub) {assert(0);}
			    
    /*! \brief Load in a problem by assuming ownership of the arguments.
	    The constraints on the rows are given by lower and upper bounds.

	For default argument values see the matching loadProblem method.

	\warning
	The arguments passed to this method will be freed using the
	C++ <code>delete</code> and <code>delete[]</code> functions. 
    */
    virtual void assignProblem (CoinPackedMatrix*& matrix,
			        double*& collb, double*& colub, double*& obj,
			        double*& rowlb, double*& rowub) {assert(0);}

    /*! \brief Load in a problem by copying the arguments.
	    The constraints on the rows are given by sense/rhs/range triplets.
	    
	If a pointer is 0 then the following values are the default:
	<ul>
          <li> <code>colub</code>: all columns have upper bound infinity
          <li> <code>collb</code>: all columns have lower bound 0 
	  <li> <code>obj</code>: all variables have 0 objective coefficient
          <li> <code>rowsen</code>: all rows are >=
          <li> <code>rowrhs</code>: all right hand sides are 0
          <li> <code>rowrng</code>: 0 for the ranged rows
        </ul>

	Note that the default values for rowsen, rowrhs, and rowrng produce the
	constraint ax >= 0.
    */
    virtual void loadProblem (const CoinPackedMatrix& matrix,
			      const double* collb, const double* colub,
			      const double* obj,
			      const char* rowsen, const double* rowrhs,   
			      const double* rowrng) {assert(0);}

    /*! \brief Load in a problem by assuming ownership of the arguments.
	    The constraints on the rows are given by sense/rhs/range triplets.
	
	For default argument values see the matching loadProblem method.

	\warning
	The arguments passed to this method will be freed using the
	C++ <code>delete</code> and <code>delete[]</code> functions. 
    */
    virtual void assignProblem (CoinPackedMatrix*& matrix,
			        double*& collb, double*& colub, double*& obj,
			        char*& rowsen, double*& rowrhs,
			        double*& rowrng) { assert(0); }

    /*! \brief Load in a problem by copying the arguments. The constraint
	    matrix is is specified with standard column-major
	    column starts / row indices / coefficients vectors. 
	    The constraints on the rows are given by lower and upper bounds.
    
      The matrix vectors must be gap-free. Note that <code>start</code> must
      have <code>numcols+1</code> entries so that the length of the last column
      can be calculated as <code>start[numcols]-start[numcols-1]</code>.

      See the previous loadProblem method using rowlb and rowub for default
      argument values.
    */
    virtual void loadProblem (const int numcols, const int numrows,
			      const CoinBigIndex * start, const int* index,
			      const double* value,
			      const double* collb, const double* colub,   
			      const double* obj,
			      const double* rowlb, const double* rowub) { assert(0); }

    /*! \brief Load in a problem by copying the arguments. The constraint
	    matrix is is specified with standard column-major
	    column starts / row indices / coefficients vectors. 
	    The constraints on the rows are given by sense/rhs/range triplets.
    
      The matrix vectors must be gap-free. Note that <code>start</code> must
      have <code>numcols+1</code> entries so that the length of the last column
      can be calculated as <code>start[numcols]-start[numcols-1]</code>.

      See the previous loadProblem method using sense/rhs/range for default
      argument values.
    */
    virtual void loadProblem (const int numcols, const int numrows,
			      const CoinBigIndex * start, const int* index,
			      const double* value,
			      const double* collb, const double* colub,   
			      const double* obj,
			      const char* rowsen, const double* rowrhs,   
  			      const double* rowrng) { assert(0); }

 


  //@}

  //---------------------------------------------------------------------------

 
 
  /*! \name OsiSimplexInterface
      \brief Simplex Interface

    Methods for an advanced interface to a simplex solver. The interface
    comprises two groups of methods. Group 1 contains methods for tableau
    access. Group 2 contains methods for dictating individual simplex pivots.
  */
  //@{

  /*! \brief Return the simplex implementation level.
  
      The return codes are:
      - 0: the simplex interface is not implemented.
      - 1: the Group 1 (tableau access) methods are implemented.
      - 2: the Group 2 (pivoting) methods are implemented

      The codes are cumulative - a solver which implements Group 2 also
      implements Group 1.
  */
  virtual int canDoSimplexInterface() const { return 0; }
  

    
  //---------------------------------------------------------------------------

  ///@name Constructors and destructors
  //@{
    /// Default Constructor
    OsiSubproblemWrapper(stochasticInput& in, int whichScenario); 
    
    virtual OsiSolverInterface * clone(bool copyData = true) const { assert(0); return 0; }
  
    /// Destructor 
    virtual ~OsiSubproblemWrapper() {}

      //@}

  //---------------------------------------------------------------------------
private:
 
	CoinPackedMatrix rMat, cMat;
	std::vector<double> collb, colub, rowlb, rowub, rhs, colsol, coldualsol;
	std::vector<double> rowActivity;
	std::vector<bool> isInteger;
	std::vector<char> rowSense;
	int nvar1, nvar2;
	int ncons2;

};


#endif
