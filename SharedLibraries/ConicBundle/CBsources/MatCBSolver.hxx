/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatCBSolver.hxx

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



#ifndef CONICBUNDLE_MATCBSOLVER_HXX
#define CONICBUNDLE_MATCBSOLVER_HXX

/**  @file MatCBSolver.hxx
    @brief Header declaring the classes ConicBundle::MatrixCBSolver, ConicBundle::MatrixFunctionOracle and ConicBundle::PrimalMatrix
    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg
*/


#include "CBSolver.hxx"
#include "matrix.hxx"

//------------------------------------------------------------


namespace ConicBundle {

/**@defgroup cxxmatrixinterface Interface to ConicBundle for the Language C++ using Matrix Classes
   @brief Solve \f$min_{y\in\mathbf{R}^m}  f_0(y) + f_1(y) + ... + f_k(y)\f$ 
    for convex functions f_i, the y-variables may be bounded or 
    box constrained. The most important steps are explained here

   Before starting, note that this interface relies on several classes
   defined in CH_Matrix_Classes, in particular on
   CH_Matrix_Classes::Matrix and CH_Matrix_Classes::Indexmatrix. The
   full functionality of ConicBundle is only available with these
   classes. If you prefer an interface without having to use these,
   please use the 
   \ref cxxinterface "Interface to ConicBundle for the Language C++".  
   We now give a short overview of the most important
   steps in using the ConicBundle solver.

   <b>Setting up the Problem, the Functions, and the Main Loop</b>

   First create a new problem/solver ConicBundle::MatrixCBSolver, let
   us call it solver for brevity. In invoking the constructor a boolean
   flag may be set to enforce the use of a minimal bundle model that
   employs only one common aggregate and one common subgradient for
   all functions, so basically "no bundle", which may be favorable if
   fast iterations and/or little memory consumption are essential.
 
   Next, set the dimension of the
   design variables/argument as well as possible box constraints on
   these by ConicBundle::MatrixCBSolver::init_problem().

   Now set up each of your functions f_i as a
   ConicBundle::MatrixFunctionOracle.  Via the routine
   ConicBundle::MatrixFunctionOracle::evaluate() you will supply, for
   a given argument, the function value and a subgradient (=the
   gradient if the function is differentiable) to the solver. The
   function evaluate() is the only function that you definitely have
   to provide, see the miniature example below.

   The function oracles have to be added to the solver 
   using the routine  ConicBundle::MatrixCBSolver::add_function().

   Once all functions are added, the optimization process can be
   started. If you know a good starting point then set it with
   ConicBundle::MatrixCBSolver::set_new_center_point() now, otherwise
   the method will pick the zero vector or, in the case of box
   constraints, the point closest to zero as starting point.

   Finally, set up a loop that calls
   ConicBundle::MatrixCBSolver::do_descent_step() until
   ConicBundle::MatrixCBSolver::termination_code() is nonzero.


   <b>Retrieving Some Solution Information</b>

   After the first call to
   ConicBundle::MatrixCBSolver::do_descent_step() you can retrieve, at
   any time, the current objective value by
   ConicBundle::MatrixCBSolver::get_objval() and the argument leading
   to this value by ConicBundle::MatrixCBSolver::get_center(). For
   some screen output, use ConicBundle::MatrixCBSolver::set_out().


   <b>Lagrangean Relaxation, Primal Approximations, and Cutting Planes</b>

   If you are optimizing the Lagrange multipliers of a Lagrangean
   relaxation, you might be interested in getting an approximation to
   your primal optimal solution. This can be done by specifying in
   each function for each (epsilon) subgradient the corresponding
   primal vectors that generate it, see the parameter primal_solutions
   in ConicBundle::MatrixFunctionOracle::evaluate() as a start. Then
   for each of your functions, you can retrieve the current primal
   approximation using
   ConicBundle::MatrixCBSolver::get_approximate_primal().

   If, in addition, you plan to improve your primal relaxation
   via cutting planes, that are strongly violated by the current
   primal approximation, you should have a look at 
   ConicBundle::MatrixCBSolver::append_variables(),
   ConicBundle::MatrixFunctionOracle::subgradient_extension() and
   ConicBundle::MatrixCBSolver::reinit_function_model(). If you want to
   get rid of primal constraints/dual variables, use
   ConicBundle::MatrixCBSolver::get_approximate_slacks()
   and ConicBundle::MatrixCBSolver::delete_variables().

   @include mat_mini_ex.cxx

*/
  //@{

  /** @brief If in Lagrangean relaxation primal solutions are in the form of a real vector
or, more generally a matrix, then an approximate primal solution can be generated by supplying primal information of this form for each epsilon subgradient within ConicBundle::MatrixFunctionOracle::evaluate(). 

      In many applications, e.g. in Lagrangean relaxation, the convex
      minimization problem arises as the dual of a convex primal  
      maximization problem. In this case one is typically interested
      in obtaining a primal approximate solution in addition to the
      dual solution. Under reasonable conditions this is possible if the
      primal solutions that give rise to the subgradients are 
      aggregated along with the subgradients within the bundle algorithm.
      If the primal data can be represented as a CH_Matrix_Classes::Matrix 
      then the user has to supply in the oracle for each sugradient the 
      corresponding primal data in a PrimalMatrix and the algorithm will do the rest. 
      Observe that a PrimalMatrix can be used exactly in the same way as
      a CH_Matrix_Classes::Matrix and that they are assignable among each other.

      The primal data has to be supplied within ConicBundle::MatrixFunctionOracle::Evaluate()
      and can be retrieved via the methods 
      ConicBundle::MatrixCBSolver::get_approximate_primal() and
      ConicBundle::MatrixCBSolver::get_center_primal()
 */
 
  class PrimalMatrix:public PrimalData, public CH_Matrix_Classes::Matrix
  {

  public:
    /// empty matrix
    PrimalMatrix(){}
  /** @brief generate a matrix of size nr x nc but WITHOUT initializing the memory
      
      If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use 
      set_init() via matrix.set_init(true) in order to avoid warnings concerning improper 
      initialization
  */ 
    PrimalMatrix(CH_Matrix_Classes::Integer nr,CH_Matrix_Classes::Integer nc):Matrix(nr,nc)   {}
    /// generate a matrix of size nr x nc initializing all elements to the value d
    PrimalMatrix(CH_Matrix_Classes::Integer r,CH_Matrix_Classes::Integer c,CH_Matrix_Classes::Real d)    :Matrix(r,c,d) {}
    /// copy constructor, *this=pm
    PrimalMatrix(const PrimalMatrix& pm):PrimalData(),Matrix(pm)    {}  
    /// copy constructor, *this=pm
    PrimalMatrix(const CH_Matrix_Classes::Matrix& pm):Matrix(pm)    {} 
    /// copy operator, *this=pm
    PrimalMatrix& operator=(const CH_Matrix_Classes::Matrix& pd)  
    { Matrix::operator=(pd); return *this; }
    
    /// produces a new PrimalMatrix that is a copy of itself; the caller has to delete the returned object at some point
    PrimalData* clone_primal_data(){return new PrimalMatrix(*this);}

    /// copy its information to *this (it must dynamic_cast to a PrimalMatrix)
    int assign_primal_data(const PrimalData& it)
    {
      const PrimalMatrix* pd=dynamic_cast<const PrimalMatrix*>(&it);
      assert(pd!=0);
      *this = *pd;
      return 0;
    }

    /// multiply *this Matrix with myfactor and add itsfactor*it (it must dynamic_cast to a PrimalMatrix)
    int aggregate_primal_data(double myfactor,double itsfactor,const PrimalData& it)
    {
      const PrimalMatrix* pd=dynamic_cast<const PrimalMatrix*>(&it);
      if (pd!=0) {
	CH_Matrix_Classes::xbpeya(*this,*pd,itsfactor,myfactor);
	return 0;
      }
      return 1;
    }
  };


  /**@brief Oracle interface (abstract class). For each of your functions, provide an instance of a derived class.  

     The oracle interface is used to describe and pass convex objective 
     functions to the ConicBundle::MatrixCBSolver. 
     The dimension of the argument vector of the function must be
     set in ConicBundle::MatrixCBSolver::init_problem() and the functions
     are then added to the solver by ConicBundle::MatrixCBSolver::add_function(). 

     If the sum of several such functions is to be minimized it
     is the task of the user to guarantee, that all dimensions
     match.

     If the function corresponds to Lagrangean relaxation of a 
     primal maximization problem one may want to generate a primal 
     approximate solution. In this case, return in the function
     ConicBundle::MatrixCBSolver::evaluate() the generating primal 
     objects for each subgradient.
     If no primal objects are returned, 
     there will be no primal aggregation.

     If primal aggregation is used then it is possible to implement
     a primal cutting plane framework. This requires the introduction
     of new (dual) variables in the design space of the function. 
     In this case the function 
     ConicBundle::MatrixCBSolver::subgradient_extension() serves the
     purpose of filling in the missing coordinates in existing
     subgradients. If this feature is not needed, the function
     may be used as is and need not be reimplemented.
  */


  class MatrixFunctionOracle: public FunctionObject 
  {
  public:
        
    /** @brief Called by the solver. Has to Return function value and 
               at least one (epsilon) subgradient and, possibly for Lagrangean 
               relaxation, some primal data.

        The evaluation method is the main interface to the bundle solver. 
        The solver calls this method to obtain for the @a current_point
        (its dimension is set in ConicBundle::MatrixCBSolver::init_problem())
        the @objective_value and (epsilon) subgradient information.
        In any call several epsilon subgradients may be returned in
        @a cut_values and @a eps_subgradients, but at least one has
        to be returend. Each subgradient describes a linear minorant
        to the convex function and is used by the bundle method to form 
        a cutting model of the convex function.
       
        The i-th linear minorant consits of @a cut_values(i) (the value of the
        cutting plane in @a current_point) and epssubgradients.col(i) 
        (the i-th column of the CH_Matrix_Classes::Matrix describes 
        the linear behavior)

        In many applications, computing the function value is an iterative
        process that approaches the true function value from below. 
        On input the code offers a bound in @a objective_value for the 
        function value, above which it is certain that the code will reject 
        the current point. If in the iterative
        process a lower bound on the function value exceeds this bound, 
        then it is sufficient to return, instead of the true function value
        and a subgradient, the current lower bound and a vector so that
        together they describe a supporting hyperplane (i.e. a linear minorant) 
        to the function at this point. 

        If the function corresponds to Lagrangean relaxation of a 
        primal maximization problem one may want to generate a primal 
        approximate solution. The parameter @a primal_data serves this
        purpose. If at each call and for each epsilon subgradient the
        corresponding generating primal object (must be derived
        from ConicBundle::PrimalData, e.g., a ConicBundle::PrimalMatrix) 
        is stored in @a primal_data, then the code automatically 
        performs the aggregation corresponding to the aggregation
        of the subgradients on the primal data objects. The
        primal approximate solution is finally delivered by
        the methods ConicBundle::MatrixCBSolver::get_approximate_primal()
        or  ConicBundle::MatrixCBSolver::get_center_primal(). All primal objects
        passed to the solver via @a primal_data must be objects allocated
        on the heap. The ownership of tese objects is transferred 
        to the solver and the solver will destroy them eventually,
        so DO NOT delete them yourself!

        If no primal aggregation is desired, simply do not touch 
        @a primal_data or clear it.  

        @param[in] current_point (const Matrix&)
           argument of the function as a column vector (e.g. the Lagrange multipliers)

        @param[in] relprec (double)
           relative precision requirement for objective values 
           that may lead to descent steps (this precision is not
           required if it is already certain that the function 
           value will be too poor)

        @param[in,out]  objective_value (double&) 
         - on input: 
           value gives the threshold for a null step; you may stop, 
           if a cutting plane yields at least this; 
         - on output: 
           return an upper bound on the true function value within
           @a relprec *(abs(objval)+1.), if there is no linear minorant 
           cutting above the threshold specified in objective_value on input. 
           Otherwise the return value should be the max of cut_values.

        @param[out]  cut_values (Matrix&) 
           store for each linear minorant (epsilon subgradient) the value 
           at the argument (i-th coordinate of the cut_values vector
           must correspond to the i-th column in eps_subgradients)

        @param[out]  eps_subgradients (Matrix&)
           store for each linear minorant the epsilon subgradient vector
           as columns (column index must correspond to
           to index in cut_values).

        @param[out]  primal_data (std::vector<PrimalData*>&)
           If the function arises from Lagrangean relaxation and a 
           primal approximation is desired then return the primal data 
           objects corresponding to the eps_subgradients.
           The control over all PrimalData objects
           pointed to by primal_data is passed to the calling routine. 
           The calling routine will delete these objectes eventually.

        @param[out]  primal_extender (PrimalExtender*&)
           if primal_data of previous calls has now to be updated
           due to changes in the primal problem -- e.g., this may happen
           in column generation -- one may return a pointer to
           PrimalExtender object on the heap. This object will be used
           by ConicBundle to update all its internally stored
           primal_data objects by calling PrimalExtender::extend on
           each of these, afterwards ConicBundle deletes primal_extender
           If this is not needed, the variable holds 0.
  
        @return  
           -  0,  all correct
           -  !=0, failure. This does not necessarily terminate the bundle method.
               Termination is forced only if no new subgradient is returned. 
     */
    virtual
    int
    evaluate
    (
     const  CH_Matrix_Classes::Matrix&          current_point,
     double relprec,
     double&                 objective_value,
     CH_Matrix_Classes::Matrix&  cut_values,   
     CH_Matrix_Classes::Matrix&  eps_subgradients,
     std::vector<PrimalData*>&    primal_data,
     PrimalExtender*& primal_extender
     )
      = 0;

    /**@brief This routine is need not be implemented unless variables 
              (constraints in Lagrangean relaxation) 
              are added on the fly. 
 
         The solver calls this routine whenever new variables have been added on
         the fly in order to extend old subgradients to the new coordinates.
         If primal data was supplied for the subgradients then  
         @a generating_primal holds a pointer to this (possibly aggregated) data,
         otherwise it is NULL. 
 
         In the presence of primal data, the new coordinates correspond to 
         the violation of the new primal constraints. These have to be 
         returned in  @a new_subgradient_values; more precisely, 
         for i=0 to @a varialb_indices.size()-1 the element 
         @a new_subgradient_values[i] has to hold the subgradient information  
         of constraint @a variable_indices[i]; 
 
         If generating_primal is NULL, then the routine can only successfully
         extend the subgradients, if the new coordinates have no influence
         on the function; then the new subgradient coordinates are all zero
         and the components of @a new_subgradient_values have to be initialized 
         to zero.

        @param[in] generating_primal  (const PrimalData*, may be NULL)
            If not Null it holds the (possibly aggregated) primal solution that 
            generated the subgradient that needs to be extendend


        @param[in] variable_indices  (const Indexmatrix&)
            for the variables with indices in @a variable_indices[i], 
            the subgradient coefficient has to be computed

        @param[out] new_subgradient_values  (Matrix &)
	     store the the subgradient coefficient of the variable with index
	     @a variable_indices[i] at @a new_subgradient_values[i] for all i.
      
         @return
           -  0 on success,
           -  1 if extension is impossible 
    */

    virtual
    int
    subgradient_extension
    (
     const PrimalData* /* generating_primal */,
     const CH_Matrix_Classes::Indexmatrix& /* variable_indices */, 
     CH_Matrix_Classes::Matrix& /* new_subgradient_values */)
    {return 1;}

       
  };

  
  class MatrixBSolver;

  /**@brief   Bundle method solver.
     
  Minimizes the sum of convex functions that are given via 
  ConicBundle::MatrixFunctionOracle interfaces, see 
  \ref cxxmatrixinterface "the text explaining the C++ interface for Matrix Classes" for
  a quick overview.  

  It provides special support for Lagrangean relaxation by generating
  primal approximate solutions if such information is provided in the
  function oracles. 

  Based on these primal approximations it is also possible to implement 
  cutting plane schemes. Routines for adding and deleting corresponding
  dual variables as well as a framework for extending subgradients in order
  not to loose the cutting model are available. 
  
  */	
  class MatrixCBSolver
  {
  private:		
    MatrixBSolver* solver;                       ///< pointer to internal solver data
    MatrixCBSolver ( const CBSolver& );           ///< not available, blocked deliberately
    MatrixCBSolver& operator= ( const CBSolver& );///< not available, blocked deliberately
  public:    
    
    /** @brief if the input parameter no_bundle is set to true then a minimal bundle variant is used consisting
	of just one aggregate and one new subgradient in each iteration (this is basically "no bundle", try
        it whenever fast iterations and/or low memory consumption seem advantageous)
    */ 
    MatrixCBSolver(bool no_bundle=false);
    ///
    ~MatrixCBSolver();
    
    //------------------------------------------------------------
    /**@name Initialization */
    //@{
       
    /** @brief Clears all data structures and problem information 
        but keeps ouptut settings and algorithmic parameter settings 
    */
    void clear();

    /** @brief Sets default values for algorithmic parameters that are not function specific (e.g., relative precision, weight and weight bounds for the augmentedproblem, etc.) 
    */
    void set_defaults();

    /** @brief Initializes the problem by setting up the design space 
        (the dimension and possible box constraints of the variables)
 
    Clears all data structures and sets the dimension @ m for a new problem.
    for solving   min_{y in R^m}  f_0(y) + f_1(y) + ...
    Box constraints may be specified for y. (The functions f_i must be added 
     by add_function()). 
  
     Lower and/or upper bounds must be speicified for all variables
     or for none of them. To specify no bounds at all, give Null
     pointers. Otherwise use ConicBundle::CB_minus_infinity for 
     unbounded below and ConicBundle::CB_plus_infinity for unbounded above.
     For NULL pointers, unbounded will be used as default for all
     variables. Specifying bounds selectively is also possible
     by set_lower_bound() or set_upper_bound().

     @param[in] dimm  (int)
         the dimension of the argument/design space/the number of Lagrange multipliers
   
     @param[in] lbounds  (const Matrix*)
         If NULL, all variables are considered unbounded below,
         otherwise lowerb[i] gives the minimum feasible value for variable y[i],
         use ConicBundle::CB_minus_infinity for unbounded below.

     @param[in] ubounds (const Matrix*)
         If NULL, all variables are considered unbounded above,
         otherwise upperb[i] gives the maximum feasible value for variable y[i],
         use ConicBundle::CB_plus_infinity for unbounded above.

     @param[in] costs (const Matrix*)
         Use this in order to specify linear costs on the variables in addition
	 to the functions (may be convenient in Lagrangean relaxation for
	 the right hand side of coupling contsraints); NULL is equivalent
	 to costs zero.

     @return 
        - 0 on success
        - != 0 otherwise

 
     */
    int init_problem(int dim, 
		     const CH_Matrix_Classes::Matrix* lbounds=0,
		     const CH_Matrix_Classes::Matrix* ubounds=0,		     
		     const CH_Matrix_Classes::Matrix* costs=0);
    

    /** @brief Adds a function, typically derived from ConicBundle::FunctionOracle; all functions added must have the same argument dimension set in init_problem().

     Besides the standard ConicBundle::MatrixFunctionOracle 
     the interface only accepts a few other prespecified derivations 
     of the class FunctionObject that come along with the CH_Matrix_Classes interface
     (e.g. for semidefinite and second order cones). Functions not derived from
     these will fail to be added and return a value !=0.

    @return 
      - 0 on success
      - != 0 otherwise
    */

    int
    add_function
    ( FunctionObject& function );
    
    /**@brief Sets lower bound for variable i,
	use ConicBundle::CB_minus_infinity for unbounded from below.

       The algorithm may have to adapt the center point aftwards.
       In this case the old function values will be marked as outdated and 
       will be recomputed at the next call to e.g. do_descent_step().
 
    @return 
      - 0 on success
      - != 0 otherwise
    */
    int set_lower_bound(int i, double lb);

    /**@brief Sets upper bound for variable i,
	use ConicBundle::CB_plus_infinity for unbounded from below.

       The algorithm may have to adapt the center point aftwards.
       In this case the old function values will be marked as outdated and 
       will be recomputed at the next call to e.g. do_descent_step().

     @return 
      - 0 on success
      - != 0 otherwise
    */
    int set_upper_bound(int i, double ub); 


    /** @brief Append new variables (always in last postions in this order).

       If 0 is feasible for the new coordinates then this is selected as
       starting value for the new coordinates; otherwise, the number
       closest to zero is used. If all new coordinates can be set to zero 
       then it is assumed that for an existing center point the 
       function values need not be recomputed (this is e.g. the case in 
       Lagrangean relaxation; if this is not correct call 
       reinit_function_model() below). 
       Otherwise the old function values will be marked as outdated and 
       will be recomputed at the next call to e.g. do_descent_step().

     @attention Be sure to update your objective functions so that they can 
       handle the new variables before you call this and any further ConicBundle 
       routines that require function evaluations. 
       Also, these operations may lead to inavailability of certain 
       other data such as subgradients and primal approximations.
 
     @param[in] n_append  (int)
       number of variables to append (always in last position in the same order)
  
     @param[in] lbounds  (const Matrix*)
        If NULL, all appended variables are considered unbounded below,
        otherwise lowerb[i] gives the minimum feasible value for variable y[i],
        use ConicBundle::CB_minus_infinity for unbounded below.

     @param[in] ubounds (const Matrix*)
        If NULL, all appended variables are considered unbounded above,
        otherwise upperb[i] gives the maximum feasible value for variable y[i],
        use ConicBundle::CB_plus_infinity for unbounded above.

     @return 
        - 0 on success
        - != 0 otherwise


    */
    int append_variables(int n_append, 			 
			 const CH_Matrix_Classes::Matrix* lbounds=0,
			 const CH_Matrix_Classes::Matrix* ubounds=0,
			 const CH_Matrix_Classes::Matrix* costs=0);

    /** @brief Deletes variables corresponding to the specified indices. 

        The indices of the remaining variables are reassigned so that they
        are consecutive again, the routine returns in @a map_to_old 
        a vector giving for each new index of these remaining variables 
        the old coordinate.

        If all of the deleted variables are zero, function values are assumed 
        to remain correct (if this is not so, call reinit_function_model() below)
        Otherwise the old function values will be marked as outdated and 
        will be recomputed at the next call to e.g. do_descent_step().

     @attention Be sure to update your objective functions so that they can 
        handle the new variables before you call any further ConicBundle 
        routines that require function evaluations.
        Also, these operations may lead to inavailability of certain 
        other data such as subgradients and primal approximations.

     @param[in] delete_indices  (const Indexmatrix&)
        the entries delete_indices[i] specify the indices of the variables 
        to be deleted
  
     @param[out] map_to_old  (Indexmatrix&) 
        after the call, element map_to_old[i] gives the old index (before the call)
        of the variable that now has index position i.
  
     @return 
        - 0 on success
        - != 0 otherwise

    */
    int delete_variables(const CH_Matrix_Classes::Indexmatrix& del_indices,CH_Matrix_Classes::Indexmatrix& map_to_old);

    /** @brief Reassigns variables to new index positions by mapping to position @a i 
        the variable that previously had index @a assign_new_from_old[i].

        Old variables, that are not mapped to any position will be deleted.
        It is allowed to generate several copies of old variables. 

        If all of the deleted variables as well as new multiple copies are zero, 
        function values are assumed to remain correct (if this is not so, 
        call reinit_function_model() below).
        Otherwise the old function values will be marked as outdated and 
        will be recomputed at the next call to e.g. do_descent_step().

     @attention Be sure to update your objective functions so that they can 
        handle the new variables before you call any further ConicBundle 
        routines that require function evaluations.
        Also, these operations may lead to inavailability of certain 
        other data such as subgradients and primal approximations.

     @param[in] assign_new_from_old  (const IVector&)
        entry assign_new_from_old[i] specifies
        the old index of the variable, that has to be copied to index position i.
    
     @return 
        - 0 on success
        - != 0 otherwise

    */
    int reassign_variables(const CH_Matrix_Classes::Indexmatrix& assign_new_from_old);

    //@}

    //------------------------------------------------------------
    /**@name Basic algorithmic routines and parameters */
    //@{
    
    /** @brief Does a descent step for the current center point.

        A descent step may consist of several function evaluations
        (null steps), that lead to no immediate progress but serve 
        for building a cutting model of the objective function 
        close to the current center point.
        A minimizer to the model is accepted as descent step
        if the function value at this point satisfies a sufficient 
        decrease criterion in comparison to the decrease predicted 
        by the model. Having found a descent step, the next center is 
        automatically shifted to this successful candidate. 
        Termination criteria may stop the process of seeking for
        a descent step, in which case the current center is kept
        and the routine termination_code() returns the termination
        code. 

        Restarting, after each descent step, the bundle method from scratch 
        with the new center as starting point does not endanger convergence.
        Therefore, a descent step is the smallest unit, after which
        user interaction can take place safely and this is the default
        choice.

	If you know what your are doing, you may also use the input
        parameter maxsteps to force the algorithm to return after
        at most maxsteps null steps. Calling do_descent_step again
        without any intermediate problem configurations will then
        simply continue the process where it stopped and convergence
        is save. During null steps one may not decrease the weight
        or delete nonzero variables of the center or the current candidate!
     
        In a Lagrangean relaxation cutting plane approach one may want
        to separate and enlarge the dimension after a certain number
        of null steps. In this case the code will try to preserve the model,
        given appropriate subgradient extension routines have been 
        provided. If the model cannot be extended, it has to be
        discarded (if subgradient extension is not successful
        this is done automatically), and the algorithm will be restarted
        from the current center point. 
	 
      @param[in] maxsteps (int)
          if maxsteps>0 the code returns after at most so many null steps 

      @return 
        - 0 on success
        - != 0 otherwise
    
    */
    int 
    do_descent_step(int maxsteps=0);
    
    
    /** @brief Returns the termination code of the bundle algorithm for the latest descent step

      For resetting all counters relevant for termination see clear_fail_counts() .
 
      @return
      -  0  :    Not terminated. 
             (Continue with the next do_descent_step())
      -  1  :    Relative precision criterion satisfied. (See set_term_relprec())
      -  2  :    Timelimit exceeded. 
             (Currently the C interface does not offer a timelimit.)
      -  4  :    Maximum number of function reevaluations exceeded. 
             (Indicates that there is a problem with one of the function 
             oracles that seems to deliver no valid upper bounds on the true 
             function value for descent steps)
      -  8  :    Maximum number of quadratic subproblem failures exceeded.
             (Indicates that the numerical limits of the inner quadratic 
             programming solver are reached, no further progress expected) 
      - 16  :    maximum number of model evaluation failures exceeded
             (Indicates that the numerical limits of the setup of the 
             subproblem are reached, no further progress expected) 
      - 32  :    maximum number of failures to increase the augmented model value exceeded
             (Indicates that the numerical limits  of the interplay between 
             subproblem and quadratic programming solver are reached, 
             no further progress expected) 
       - 64  :   maximum number of oracle calls (function evaluations) exceeded,
                 see set_eval_limit()
       - 128  :   maximum number of oracle failures exceeded. 
             This refers to function evaluations that terminate with insufficient 
	     precision but still provide a new approximate subgradient. A failure typically
             indicates numerical difficulties with the precision requirements.  
             (Currently the interface does not allow to manipulate the limit, it is set to 10)


    */
    int
    termination_code() 
      const; 
    
    /** @brief Outputs a text version of termination code, see termination_code(). 

      @return 
        - 0 on success
        - != 0 otherwise

    */
    std::ostream&
    print_termination_code(std::ostream& out); 
    
    /** @brief Returns the objective value resulting from last descent 
        step (initially undefined). If no problem modification routines 
        were called since then, it is the objective value at the point 
        returned by get_center().
    */
    double 
    get_objval()
      const;
    
    /** @brief Returns the next center point that was produced by the latest call  
	to do_descent_step (in some problem modification routines the
	center point may be updated immediately, in others the center point 
	will be corrected automatically directly before starting 
	the next descent step and its values may be infeasible till then).

      @return 
        - 0 on success
        - != 0 otherwise
    */
    int 
    get_center
    ( CH_Matrix_Classes::Matrix& center )
      const;
    
    
    /** @brief Returns Euclidean norm of the latest aggregate subgradient.
     */
    double 
    get_sgnorm()
      const; 

    /** @brief Returns the latest aggregate subgradient.

    @return 
      - 0 on success
      - != 0 otherwise

    */
    int 
    get_subgradient
    ( CH_Matrix_Classes::Matrix& subgradient )
      const;

    /** @brief Returns the cutting model value resulting from last call to 
        do_descent_step() (initially undefined). 
    */
    double 
    get_cutval()
      const;
    
    /** @brief Returns the objective value computed in the last step of do_descent_step(), 
        independent of whether this was a descent step or a null step (initially undefined). 

        If no problem modification routines were called since then, it is the 
        objective value at the point returned by get_candidate(). If this 
        last evaluation led to a descent step, then it is the same value as
        in get_objval().
    */
    double 
    get_candidate_value()
      const;
    
    /** @brief Returns the last point, the "candidate", at which the function 
        was evaluated in do_descent_step(). 

        If this evaluation lead to a descent step, it is the same point as 
	in get_center().

      @return 
        - 0 on success
        - != 0 otherwise
    */
    int 
    get_candidate
    ( CH_Matrix_Classes::Matrix& center )
      const;
    
    
    //@}
    
    //------------------------------------------------------------
    /**@name Advanced algorithmic routines and parameters */
    //@{

    /** @brief Sets the relative precision requirements for successful termination
               (default 1e-5).

     @param[in] term_relprec (double)
       The algorithm stops with termination code 1, if predicted progress for 
       the next step is less than term_relprec times
       absolute function value plus one. 

    @return 
      - 0 on success
      - != 0 otherwise

    */
    int 
    set_term_relprec
    ( const double term_relprec );
    
    /** @brief Set the starting point/center that will be used in the 
        next call to  do_descent_step(). Each call
        to this routine causes an immediate evaluation of all oracles. 

     @return  
        - 0 on success 
        - != 0 otherwise
    */
    int
    set_new_center_point
    ( const CH_Matrix_Classes::Matrix& center_point );

    
    /** @brief Returns the return value of the latest evaluation call
        to this @a function.
    */
     int 
    get_function_status
    ( const FunctionObject& function )
      const;

     /** @brief Returns the multipliers for the box constraints on the design variables;
        in Lagrangean relaxation they may be interpreted as primal slacks
	for inequality constraints.
     @return 
        - 0 on success
        - != 0 otherwise
    */
    int 
    get_approximate_slacks(CH_Matrix_Classes::Matrix &)
      const;
    
    /** @brief returns the current approximate primal solution corresponding
        to the aggregate subgradient of the specified @a function. 
        
        PrimalData solutions must have been supplied in all previous 
        calls to evaluate; In this case it returns the current approximate 
        primal solution aggregated alongside with the aggregate subgradient.  
        A primal solution may not be available after addition of constraints, 
        if extension of the aggregate subgradient to the new coordinates failed.

     @return 
        - 0 on success
        - != 0 otherwise
    */
    int
    get_approximate_primal
    ( const FunctionObject& function,
      PrimalData&              primal )
      const; 
    
    /** @brief Returns the primal solution corresponding to the best epsilon
	subgradient returned in the evaluation of the specified @a function 
        at the current center point 
        
        PrimalData solutions must have been supplied in all previous 
        calls to evaluate; It may not be available or may correspond to an 
        aggregate primal after addition or deletion of design variables/primal 
        constraints.
  
     @return 
        - 0 on success
        - != 0 otherwise
    */
    int
    get_center_primal
    ( const FunctionObject& function,
      PrimalData&              primal )
      const; 


    /** @brief Returns the primal solution returned by the last
	evaluation of the specified @a function 
        in the point get_candidate().
         
        It will only be available if also supplied by the @a function 
       
     @return 
        - 0 on success
        - != 0 otherwise
    */
    int
    get_candidate_primal
    ( const FunctionObject& function,
      PrimalData&              primal )
      const; 


    /** @brief Sets the maximum number of subgradients used in forming the
        cutting model of the specified @a function  

        Quite often a very small model, e.g., 2, yields very fast iterations
        and good progress in time (sometimes at the cost of more evaluations).
        By limited numerical experience, a significant reduction in the number of 
        evaluations can  only be expected if the bundle is large enough to 
        wrap the function rather tightly. Quite frequently, unfortunately,
        this entails that solving the quadratic subproblems
        is more expensive than function evaluation. 

        The meaning of this routine may differ from standard for 
	predefined special functions with special bundle types.

     @param[in] function (const FunctionObject&)
        the function added in add_function()

     @param[in] bundlesize (int)
         maximum number of subgradients to be used in forming the cutting model 

     @return 
        - 0 on success
        - != 0 otherwise

    */
    int 
    set_max_bundlesize
    (  const FunctionObject& function, int max_bundlesize );

    /** @brief Sets the maximum number of new subgradients to be used in the
	next bundle update of the cutting modle for the specified @function.

        The meaning of this routine may differ from standard for 
	predefined special functions with special bundle types.

     @param[in] function (const FunctionObject&)
        the function added in add_function()

     @param[in] max_new_subgradients (int)
       maximum number of new epsilon subgradients to be used in bundle updates 

     @return 
       - 0 on success
       - != 0 otherwise

    */
    int 
    set_max_new_subgradients
    (  const FunctionObject& function, int max_new_subgradients );

    /**@brief Sets the maximum bundlesize and the maximum number of new subgradients
        added in a bundle update of the cutting model for the specified @a function.
        The meaning of this routine may differ from standard for 
	predefined special functions with special bundle types.

     @param[in] function (const FunctionObject&)
        the function added in add_function()

     @param[in] params (const BundleParameters&)
       some update parameters for the cutting model, see e.g. ConicBundle::BundleParameters

     @return 
       - 0 on success
       - != 0 otherwise

    */
    int 
    set_bundle_parameters
    (  const FunctionObject& function,
       const BundleParameters& params );

     /** @brief Retrieves current bundle parameters (not the actual size in use!)
        as set for the cutting model of the specified @a function.

	This may differ for predefined special 
	functions with derived BundleParameter classes. 

     @param[in] function (const FunctionObject&)
        the function added in add_function()

     @param[out] params (BundleParameters&)

     @return 
       - 0 on success
       - != 0 otherwise

    */
    int 
    get_bundle_parameters
    ( const FunctionObject& function, BundleParameters& params )
      const;

    /** @brief Returns the current bundle values: the current bundle_size and the
        number of subgradients added in the latest update of the cutting
        model of the specified @a function.
  
	This may differ for predefined special 
	functions with derived BundleParameter classes. 

     @param[in] function (const FunctionObject&)
        the function added in add_function()

     @param[out] params (BundleParameters&)

     @return 
       - 0 on success
       - != 0 otherwise

    */ 
    int 
    get_bundle_values
    ( const FunctionObject& function, BundleParameters& params )
      const;

    /** @brief Clears cutting model, subgradients and stored function values 
       for the specified @a function
 
        This has to be called whenever the specified function was modified
        so that the old subgradients and/or primal generators are no longer
        valid. 

      @param[in] function (const FunctionObject&)
        the function added in add_function()

      @return 
        - 0 on success
        - != 0 otherwise

    */
    int reinit_function_model
    ( const FunctionObject& function );
    
    /** @brief Clears the aggregate parts of the cutting model of this @a function.
 
        This has to be called whenever the specified function was modified
        so that the old aggregate subgradients and/or primal generators are 
        no longer valid. 

      @param[in] function (const FunctionObject&)
        the function added in add_function()

      @return 
        - 0 on success
        - != 0 otherwise

    */
    int clear_aggregates
    ( const FunctionObject& function );
    

    /** @brief Asks @a function to call @a primal_extender for each of its primal objects (see
     also FunctionOracle::evaluate() )

     If the function is the Lagrangian dual of a primal problem and primal_data 
     returned previous calls to the oracle has now to be updated due to changes 
     in the primal problem -- e.g., this may happen in column generation -- the 
     call causes updates of all internally stored primal_data objects by calling 
     PrimalExtender::extend on each of these.

      @param[in] function (const FunctionObject&)
        the function added in add_function()

      @param[in] primal_extender (PrimalExtender&)
        the object holding the extension function for primal_data

      @return 
        - 0 on success
        - 1 if for this function it is not possible to use a primal_extender
        - 2 if the primal_extender would be applicable but there is no primal_data
    */
    int call_primal_extender
    (const FunctionObject& function,
     PrimalExtender& primal_extender);


    /** @brief Returns the current weight for the quadratic term in the augmented subproblem
    (may be interpreted as 1./step_size or 1./trustregion-radius).
    */
    double 
    get_last_weight()
      const; 

    /** @brief Sets the  weight (>0) to be used in the quadratic term 
        of the next augmented subproblem 
        (may be interpreted as 1./step_size or 1./trustregion-radius).

        Independent of whether the weight violates current min- and max-bounds 
        set in set_min_weight() and set_max_weight(), the next model will 
        be computed for this value. Thereafter, however, it will be updated as 
        usual; in particular, it may be truncated by min and max bounds 
        immediately after the first subproblem. 

        In order to guarantee a constant weight (e.g. 1 is frequently a reasonable 
        choice if the automatic default heuristic performs poorly), set the min and max
         bounds to the same value, too.

      @param[in] weight (double)

      @return 
        - 0 on success
        - != 0 otherwise
    */
    int 
    set_next_weight
    ( const double weight );
    
    /** @brief Sets a lower bound on the  weight for the quadratic term of the 
        augmented subproblem.

        Nonpositive values indicate no bound. 
        The new value shows its effect only at first dynamic change of 
        the weight.

     @param[in] min_weight (double)

     @return 
        - 0 on success
        - != 0 otherwise

    */
    int 
    set_min_weight
    ( const double min_weight ); 

    /** @brief Sets an upper bound on the  weight for the quadratic term of the 
        augmented subproblem.

        Nonpositive values indicate no bound. 
        The new value shows its effect only at first dynamic change of 
        the weight.

      @param[in] max_weight (double)

      @return 
        - 0 on success
        - != 0 otherwise
    */
    int 
    set_max_weight
    ( const double max_weight );
    
    /** @brief Adjusts on all conic functions the penalty parameter for 
	conic violations to twice the trace of the primal approximation.

        This routine is only needed for conic function objects such
        as the nonnegative cone, the second order cone and 
        the semidefinite cone if no good upper bound on the trace of
        feasible points is known and has to be determined automatically.

        If after some time, the trace values settle, the upper bounds
        on the trace may be way to high and can then be reset with this
        call.

      @return 
        - 0 on success
        - != 0 otherwise
    */
    int 
    adjust_multiplier
    ( void );
    

    /** @brief Use a scaling heuristic or switch off scaling alltogether.

      @return 
        - 0 on success
        - != 0 otherwise
    */
    int
    set_scaling
    ( bool do_scaling );

    /** @brief  user defined diagonal scaling,
        values greater than 1 allow more movement for this variable,
        values smaller than 1 allow less movement. 

        It is the users responsibility to guarantee that the
        scaling vector fits in dimension to the runnig problem
        data, in particular if routines such as append_variables(),
        delete_variables(), and reassign_variables() are used.

      @param[in] scale (const Matrix&)

      @return 
        - 0 on success
        - != 0 otherwise
    */
    int
    set_scaling
    ( const CH_Matrix_Classes::Matrix& scale );

    /** @brief If set to true (the default is false),  
        some variables will be fixed automatically to the center 
        value if their bounds are strongly active (i.e., the
        corresponding multipliers are big). 

        The coordinates to be fixed are redetermined in each
        call following a descent step or a change of the function. 
        An indicator vector of the variables fixed in the last call 
        can be obtained via the routine get_active_bounds_indicator().

        Setting this value to true might improve the performance
        of the algorithm in some instances but there is no 
        convergence theory. It might be particularly helpful
        within Lagrangian relaxation if a primal cutting plane
        approach is used and non-tight inequalities should be
        eliminated quickly (fixing then indicates large primal 
        slack values). 

      @param[in] allow_fixing (bool)

      @return 
        - 0 on success
        - != 0 otherwise
    */ 
    void
    set_active_bounds_fixing
    ( bool allow_fixing );

    /** @brief clears all fail counts on numerical function oder model failures,
	may be useful if this caused premature termination. 

      @return 
        - 0 on success
        - != 0 otherwise
    */
    void 
    clear_fail_counts
    (void);
    
    /** @brief Sets an upper bound on the number of calls to the oracle (use negative numbers for no limit). 

        If this number is reached, the algorithm will terminate 
	independently of whether the last step was a descent or 
	a null step. A negative number will be interepreted as 
	no limit.
 
      @param[in] eval_limit (Integer)
      
      @return 
        - 0 on success
        - != 0 otherwise
    */
    void 
    set_eval_limit
    (int eval_limit);

    /** @brief Set an upper bound on the number of inner updates for the
        cutting model with primal slacks within one null step (use negative numbers for no limit). 

        A negative number will be interepreted as no limit, i.e., 
        the updates will be done till a certain precision of the 
        cutting model is achieved.
 
      @param[in] update_limit (Integer)
      
      @return 
        - 0 on success
        - != 0 otherwise
    */
    void 
    set_inner_update_limit
    (int update_limit);
    
    //@}
    
    //------------------------------------------------------------
    /**@name Look up basic paramaters (dimension, number of functions, ...)*/
    //@{
    /** @brief Returns the current dimension of the design space/argument
           or -1 if no dimension is set. 
    */
    int get_dim() const;

    /** @brief Returns the current number of functions in the problem. 
    */
    int get_n_functions() const;

    /** @brief Returns the number of function evaluations
    */
    int get_n_oracle_calls() const;

    /** @brief Returns the number of function descent setps
    */
    int get_n_descent_steps() const;

    /** @brief Returns the number of inner iterations of the bundle method
    */
    int get_n_inner_iterations() const;

    /** @brief Returns the number of inner multiplier updates for the box constraints
    */
    int get_n_inner_updates() const;

    /** @brief Returns the vector of lower bounds
    */
    const CH_Matrix_Classes::Matrix& get_lbounds() const;

    /** @brief Returns the vector of upper bounds
    */
    const CH_Matrix_Classes::Matrix& get_ubounds() const;
    
    /** @brief Returns the indicator vector of variables temporarily fixed to 
	the center value due to significantly positive multipliers
        for the box constraints.

        Such a fixing indicates that the corresponding 
        variables would like to stay at their bounds.
	If no variables were fixed, the dimension of the vector is zero.
     */
    const CH_Matrix_Classes::Indexmatrix& get_active_bounds_indicator() const;

    //@}

    //------------------------------------------------------------
    /**@name Output */
    //@{

    /** @brief Specifies the output level (out==NULL: no output at all, 
           out!=NULL and level=0: errors and warnings, 
           level>0 increasingly detailed information)

     @param[in] out  (std::ostream*) 
       direct all output to (*out). If out==NULL, there will be no output at all.

     @param[in] print_level (int)

     Output levels for print_level: 
      -  0 ... no output except for errors and warnings 
      -  1 ... line summary after each descent step  
      - >1 ... undocumented and increasingly detailed log information. 
             These higher levels should only be used if requested 
             for debugging purposes. 

      Example for level 1:

\verbatim
00:00:00.00 endit  1   1   1   563.    563.  39041.188  39043.162
00:00:00.00 endit  2   2   2   563.    559.  38488.165  38490.200
00:00:00.00 endit  3   3   3   56.3    555.  33014.533  33211.856
00:00:00.00 endit  4   4   4   5.63    517. -14306.459  2738.0343
00:00:00.00 endit  5   5   5   4.04    148. -2692.1131  2.2150883
00:00:00.00 endit  6   6   6   4.01    1.29  1.7908952  2.0000581
00:00:00.00 endit  7   7   7   3.95  0.0213  1.9999387  2.0000000
00:00:00.00 _endit  8   8   8   3.95 2.94e-05  2.0000000  2.0000000

Column 1      2     3   4   5    6       7       8          9
\endverbatim
      - Column 1: computation time in hh:mm:ss.dd,
      - Column 2: "endit" is convenient for grep and stands for "end of iteration".
         Iterations with termination_code()!=0 are marked with "_endit".
      - Column 3: number of descent steps (= calls to do_descent_step())
      - Column 4: number of descent and null steps. Up to initialization calls
         and reevaluations, this is the number of evaluation calls 
         to the function oracles from within the bundle method. 
         In the example all calls led to descent steps.
      - Column 5: number of innermost iterations. It differs from column 5 only in the
          case of variables with bounds in which case it gives the number of updates
          of the multipliers for the bounds (or primal slacks in Lagrangean 
          relaxation). Exceedingly high numbers in this column indicate that
          some variables are constantly at their bounds and it might be 
          possible to improve convergence by deleting them (i.e. set them
          as constants to this bound and remove the variable). 
      - Column 6: the weight of the quadratic term in the augmented problem.
      - Column 7: the norm of the aggregate subgradient. If it is small,
          say below 0.1, then mostly this is good indication that the
          objective value is close to optimal. 
      - Column 8: the value of the cutting model in the last candidate point. It
          is always a lower bound on the true function value in this point
      - Column 9: the objective value in the latest point that led to a descent 
          step, i.e., the point returend by get_center(). Whenever 
          termination_code() returns 0 this is also the objective
          value of the latest evaluation call to the function oracles and
          the value in the center point of the next iteration. 
     */
    void 
    set_out
    (std::ostream* out=0,int print_level=1);

    /// print a one line summary of important evaluation data
    std::ostream& print_line_summary(std::ostream& out) const;

    //@}
    
  };

  //@}
}

#endif



