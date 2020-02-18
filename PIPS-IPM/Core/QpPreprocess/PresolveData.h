/*
 * PresolveData.h
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_
#define PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_

#include "sData.h"
#include "StochPostsolver.h"
#include "StochVectorHandle.h"
#include "SimpleVectorHandle.h"
#include "SparseStorageDynamic.h"
#include "SystemType.h"

#include <algorithm>
#include <list>
#include <limits>
#include <queue>

class PresolveData
{
private:
      sData* presProb;

      StochPostsolver* const postsolver;

      const int length_array_outdated_indicators;
      bool* array_outdated_indicators;
      bool& outdated_lhsrhs;
      bool& outdated_nnzs;
      bool& outdated_linking_var_bounds;
      bool& outdated_activities;
      bool& outdated_obj_vector;
      bool& postsolve_linking_row_propagation_needed;


      /* counter to indicate how many linking row bounds got changed locally and thus need activity recomputation */
      int linking_rows_need_act_computation;

      /* number of non-zero elements of each row / column */
      StochVectorBaseHandle<int> nnzs_row_A;
      StochVectorBaseHandle<int> nnzs_row_C;
      StochVectorBaseHandle<int> nnzs_col;

      /* size of non-zero changes array = #linking rows A + #linking rows C + # linking variables */
      int length_array_nnz_chgs;
      int* array_nnz_chgs;
      SimpleVectorBaseHandle<int> nnzs_row_A_chgs;
      SimpleVectorBaseHandle<int> nnzs_row_C_chgs;
      SimpleVectorBaseHandle<int> nnzs_col_chgs;

      /* In the constructor all unbounded entries will be counted.
       * Unbounded entries mean variables with non-zero multiplier that are unbounded in either upper or lower direction.
       * Activities will be computed once the amount of unbounded variables in upper or lower direction falls below 2 so
       * that bound strengthening becomes possible.
       */
      /* StochVecs for upper and lower activities and unbounded entries */
      StochVectorHandle actmax_eq_part;
      StochVectorHandle actmin_eq_part;

      StochVectorBaseHandle<int> actmax_eq_ubndd;
      StochVectorBaseHandle<int> actmin_eq_ubndd;

      StochVectorHandle actmax_ineq_part;
      StochVectorHandle actmin_ineq_part;

      StochVectorBaseHandle<int> actmax_ineq_ubndd;
      StochVectorBaseHandle<int> actmin_ineq_ubndd;

      /// changes in boundedness and activities of linking rows get stored and synchronized
      int lenght_array_act_chgs;
      double* array_act_chgs;
      SimpleVectorHandle actmax_eq_chgs;
      SimpleVectorHandle actmin_eq_chgs;
      SimpleVectorHandle actmax_ineq_chgs;
      SimpleVectorHandle actmin_ineq_chgs;

      int* array_act_unbounded_chgs;
      SimpleVectorBaseHandle<int> actmax_eq_ubndd_chgs;
      SimpleVectorBaseHandle<int> actmin_eq_ubndd_chgs;
      SimpleVectorBaseHandle<int> actmax_ineq_ubndd_chgs;
      SimpleVectorBaseHandle<int> actmin_ineq_ubndd_chgs;

      /* handling changes in bounds */
      int lenght_array_bound_chgs;
      double* array_bound_chgs;
      SimpleVectorHandle bound_chgs_A;
      SimpleVectorHandle bound_chgs_C;

      /* storing so far found singleton rows and columns */
      std::queue<INDEX> singleton_rows;
      std::queue<INDEX> singleton_cols;

      /* SimpleVectors indicating which linking rows propagated bounds and thus need to be stored */
//      SimpleVectorBaseHandle<int> eq_linking_row_propagated_bound;
//      SimpleVectorBaseHandle<int> ineq_linking_row_propagated_bound;
//
//      SimpleVectorBaseHandle<int> linking_var_bound_implied_by_linking_row;

      const int my_rank;
      const bool distributed;

      // number of children
      const int nChildren;

      // objective offset created by presolving
      double objOffset;
      double obj_offset_chgs;
      SimpleVectorHandle objective_vec_chgs;

      // store free variables which bounds are only implied by bound tightening to remove bounds later again
      StochVectorBaseHandle<int> lower_bound_implied_by_system;
      StochVectorBaseHandle<int> lower_bound_implied_by_row;
      StochVectorBaseHandle<int> lower_bound_implied_by_node;

      // TODO a vector of INDEX would be nicer
      StochVectorBaseHandle<int> upper_bound_implied_by_system;
      StochVectorBaseHandle<int> upper_bound_implied_by_row;
      StochVectorBaseHandle<int> upper_bound_implied_by_node;

      // necessary ?
      int elements_deleted;
      int elements_deleted_transposed;

public :

      PresolveData(const sData* sorigprob, StochPostsolver* postsolver);
      ~PresolveData();

      const sData& getPresProb() const { return *presProb; };

      double getObjOffset() const { return objOffset; };
      int getNChildren() const { return nChildren; };

      void getRowActivities( const INDEX& row, double& max_act, double& min_act, int& max_ubndd, int& min_ubndd) const;
      double getRowCoeff( const INDEX& row, const INDEX& col ) const;

      const StochVectorBase<int>& getNnzsRow(SystemType system_type) const { return (system_type == EQUALITY_SYSTEM) ? *nnzs_row_A : *nnzs_row_C; }
      const StochVectorBase<int>& getNnzsRowA() const { return *nnzs_row_A; }; // todo maybe this is a problem - these counters might not be up to date
      const StochVectorBase<int>& getNnzsRowC() const { return *nnzs_row_C; };
      const StochVectorBase<int>& getNnzsCol() const { return *nnzs_col; };

      int getNnzsRow(const INDEX& row) const;
      int getNnzsCol(const INDEX& col) const;

      std::queue<INDEX>& getSingletonRows() { return singleton_rows; };
      std::queue<INDEX>& getSingletonCols() { return singleton_cols; };

      sData* finalize();

      /* reset originally free variables' bounds to +- inf iff their current bounds are still implied by the problem */
      void resetOriginallyFreeVarsBounds( const sData& orig_prob );

      /* whether or not there is currently changes buffered that need synchronization among all procs */
      bool reductionsEmpty();

      /* checks activities, non-zeros and root node */
      bool presDataInSync() const;

      /// synchronizing the problem over all mpi processes if necessary
      void allreduceLinkingVarBounds();
      void allreduceAndApplyLinkingRowActivities();
      void allreduceAndApplyNnzChanges();
      void allreduceAndApplyBoundChanges();
      void allreduceAndApplyObjVecChanges();
      void allreduceObjOffset();

      /// postsolve sync events that need to be set
      void putLinkingVarsSyncEvent();

      bool wasColumnRemoved( const INDEX& col ) const;
      bool wasRowRemoved( const INDEX& row ) const;

      /// interface methods called from the presolvers when they detect a possible modification
      // todo make bool and give feedback or even better - return some enum maybe?
      void fixColumn( const INDEX& col, double value);
      void fixEmptyColumn( const INDEX& col, double val);

      void removeSingletonRow(const INDEX& row, const INDEX& col, double xlow_new, double xupp_new, double coeff);

      void syncPostsolveOfBoundsPropagatedByLinkingRows();
      bool rowPropagatedBoundsNonTight( const INDEX& row, const INDEX& col, double xlow_new, double xupp_new, double coeff_var);
      bool rowPropagatedBounds( const INDEX& row, const INDEX& col, double ubx, double lbx);

      void substituteVariableParallelRows(const INDEX& row1, const INDEX& row2, const INDEX& col1, const INDEX& col2, double scalar, double translation);
      void removeRedundantRow( const INDEX& row );
      void removeParallelRow( const INDEX& row );
      void fixColumnInequalitySingleton( const INDEX& col, double value, double coeff );
      void removeImpliedFreeColumnSingletonEqualityRow( const INDEX& row, const INDEX& col);

      void adaptObjectiveSubstitutedRow( const INDEX& row, const INDEX& col );

      // todo : hackish functions not properly working with presolve/postsolve
      void tightenRowBoundsParallelRow(SystemType system_type, int node, int row, double lhs, double rhs, bool linking);
      void tightenVarBoundsParallelRow(SystemType system_type, int node, int row, int col, bool linking);

      /* call whenever a single entry has been deleted from the matrix */
      void deleteEntry(SystemType system_type, int node, BlockType block_type, int row, int& col_idx, int& row_end);
      void updateTransposedSubmatrix( SystemType system_type, int node, BlockType block_type, std::vector<std::pair<int, int> >& elements);

      /* methods for verifying state of presData or querying the problem */
      bool verifyNnzcounters() const;
      bool verifyActivities() const;
      bool elementsDeletedInTransposed() const { return elements_deleted == elements_deleted_transposed; };

      bool nodeIsDummy(int node) const;
      bool hasLinking(SystemType system_type) const;

      bool varBoundImpliedFreeBy( bool upper, const INDEX& col, const INDEX& row);
private:
      INDEX getRowMarkedAsImplyingColumnBound(const INDEX& col, bool upper_bound);
      void markRowAsImplyingColumnBound(const INDEX& col, const INDEX& row, bool upper_bound);

      void varboundImpliedFreeFullCheck(bool& upper_implied, bool& lower_implied, const INDEX& col, const INDEX& row) const;

      /// methods for printing debug information
      // initialize row and column nnz counter
      void initNnzCounter(StochVectorBase<int>& nnzs_row_A, StochVectorBase<int>& nnzs_row_C, StochVectorBase<int>& nnzs_col) const;
      void initSingletons();

      void setUndefinedVarboundsTo(double value);

      void addActivityOfBlock( const SparseStorageDynamic& matrix, SimpleVector& min_partact, 
            SimpleVectorBase<int>& unbounded_min, SimpleVector& max_partact,
            SimpleVectorBase<int>& unbounded_max, const SimpleVector& xlow, 
            const SimpleVector& ixlow, const SimpleVector& xupp, 
            const SimpleVector& ixupp) const ;

      long resetOriginallyFreeVarsBounds(const SimpleVector& ixlow_orig, const SimpleVector& ixupp_orig, int node);

      void adjustMatrixRhsLhsBy(SystemType system_type, int node, bool linking, int row_index, double value);
      /// methods for modifying the problem
      void adjustRowActivityFromDeletion(SystemType system_type, int node, BlockType block_type, int row, int col, double coeff);
      /// set bounds if new bound is better than old bound
      bool updateUpperBoundVariable( const INDEX& col, double xupp_new)
      { return updateBoundsVariable( col, INF_NEG_PRES, xupp_new ); };
      bool updateLowerBoundVariable( const INDEX& col, double xlow_new)
      { return updateBoundsVariable( col, xlow_new, INF_POS_PRES); };

      bool updateBoundsVariable( const INDEX& col, double xlow_new, double xupp_new );
      void updateRowActivities( const INDEX& col, double xlow_new, double xupp_new, double xlow_old, double xupp_old);

      void updateRowActivitiesBlock(SystemType system_type, int node, BlockType block_type, int col,
      		 double xlow_new, double xupp_new, double xlow_old, double xupp_old);

      void updateRowActivitiesBlock(SystemType system_type, int node, BlockType block_type, int col, double bound,
            double old_bound, bool upper);

      /* compute and update activities */
      void recomputeActivities() { recomputeActivities(false); }

      /* computes all row activities and number of unbounded variables per row
       * If there is more than one unbounded variable in the min/max activity of a row
       * +/-infinity() is stored. Else the actual partial activity is computed and stored.
       * For rows with one unbounded variable we store the partial activity without that
       * one variable, for rows with zero unbounded vars the stored activity is the actual
       * activity of that row.
       */
      void recomputeActivities(bool linking_only);

      void recomputeActivities(bool linkinig_only, StochVector& actmax_eq_part, StochVector& actmin_eq_part, StochVectorBase<int>& actmax_eq_ubndd,
         StochVectorBase<int>& actmin_eq_ubndd, StochVector& actmax_ineq_part, StochVector& actmin_ineq_part, StochVectorBase<int>& actmax_ineq_ubndd,
         StochVectorBase<int>& actmin_ineq_ubndd) const;

      double computeLocalLinkingRowMinOrMaxActivity(SystemType system_type, int row, bool upper) const;
      void computeRowMinOrMaxActivity(SystemType system_type, int node, bool linking, int row, bool upper);

      void removeColumn(const INDEX& col, double fixation);
      void removeColumnFromMatrix(SystemType system_type, int node, BlockType block_type, int col, double fixation);
      void removeRow( const INDEX& row );
      void removeRowFromMatrix(SystemType system_type, int node, BlockType block_type, int row);

      void reduceNnzCounterRow(SystemType system_type, int node, bool linking, int row_index, int amount);
      void reduceNnzCounterColumn(int node, BlockType block_type, int col_index, int amount);

      /// methods for querying the problem in order to get certain structures etc. todo: move?
      StochGenMatrix& getSystemMatrix(SystemType system_type) const;
      SparseGenMatrix* getSparseGenMatrix(SystemType system_type, int node, BlockType block_type) const;

      void checkBoundsInfeasible(const INDEX& col, double xlow_new, double xupp_new) const;
public:
      void writeRowLocalToStreamDense(std::ostream& out, const INDEX& row) const;
private:
      void writeMatrixRowToStreamDense(std::ostream& out, const SparseGenMatrix& mat, int node, int row, const SimpleVector& ixupp, const SimpleVector& xupp,
            const SimpleVector& ixlow, const SimpleVector& xlow) const;
      void printVarBoundStatistics(std::ostream& out) const;

      StochVectorHandle getRowAsStochVector(SystemType system_type, int node, int row, bool linking_row) const;
      StochVectorHandle getColAsStochVector(SystemType system_type, int node, int col) const;


};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_ */

