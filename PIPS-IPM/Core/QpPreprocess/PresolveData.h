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

struct sCOLINDEX
{
      sCOLINDEX(int node, int index) : node(node), index(index) {};
      int node;
      int index;
};

/* here we use node == -2 for the linking block in the root node */
struct sROWINDEX
{
      sROWINDEX(int node, int index, SystemType system_type) :
         node(node), index(index), system_type(system_type){};
      int node;
      int index;
      SystemType system_type;
};

class PresolveData
{
private:
      sData* presProb;



      StochPostsolver* const postsolver;

      bool outdated_activities;
      bool outdated_lhsrhs;
      bool outdated_nnzs;
      bool outdated_linking_var_bounds;
// todo make all private - forbid changing of these from external sources - only return const references / const pointers
      /* number of non-zero elements of each row / column */
      StochVectorHandle nnzs_row_A;
      StochVectorHandle nnzs_row_C;
      StochVectorHandle nnzs_col;

      /* number of removed elements of linking rows/cols stored for update */
      int length_array_nnz_chgs;
      double* array_nnz_chgs;
      SimpleVectorHandle nnzs_row_A_chgs;
      SimpleVectorHandle nnzs_row_C_chgs;
      SimpleVectorHandle nnzs_col_chgs;

      /* activities of rows */
      StochVectorHandle actmax_eq;
      StochVectorHandle actmin_eq;

      StochVectorHandle actmax_ineq;
      StochVectorHandle actmin_ineq;

      /* changes in row activities */
      int lenght_array_act_chgs;
      double* array_act_chgs;
      SimpleVectorHandle actmax_eq_chgs;
      SimpleVectorHandle actmin_eq_chgs;
      SimpleVectorHandle actmax_ineq_chgs;
      SimpleVectorHandle actmin_ineq_chgs;

      /* handling changes in bounds */
      int lenght_array_bound_chgs;
      double* array_bound_chgs;
      SimpleVectorHandle bound_chgs_A;
      SimpleVectorHandle bound_chgs_C;


      /* storing so far found singleton rows and columns */
      std::list<sROWINDEX> singleton_rows;
      std::list<sCOLINDEX> singleton_cols;

      int my_rank;
      bool distributed;

      // number of children
      int nChildren;

      // objective offset created by presolving
      double objOffset;
      double obj_offset_chgs;

      int elements_deleted;
      int elements_deleted_transposed;

public :
      const sData& getPresProb() const { return *presProb; };
      const StochVector& getActMaxEq() const { return *actmax_eq; };
      const StochVector& getActMaxIneq() const { return *actmax_ineq; };
      const StochVector& getActMinEq() const { return *actmin_eq; };
      const StochVector& getActMinIneq() const { return *actmin_ineq; };
      const StochVector& getNnzsRowA() const { return *nnzs_row_A; }; // todo maybe this is a problem - these counters might not be up to date
      const StochVector& getNnzsRowC() const { return *nnzs_row_C; };
      const StochVector& getNnzsCol() const { return *nnzs_col; };

      /* methods for initializing the object */
   private:
      // initialize row and column nnz counter
      void initNnzCounter();
      void initSingletons();
      void setUndefinedVarboundsTo(double value);

   public:
      PresolveData(const sData* sorigprob, StochPostsolver* postsolver);
      ~PresolveData();

      sData* finalize();

      /* compute and update activities */
      void recomputeActivities() { recomputeActivities(false); }
      void recomputeActivities(bool linking_only);
   private:
      void addActivityOfBlock( const SparseStorageDynamic& matrix, SimpleVector& min_activities, SimpleVector& max_activities,
            const SimpleVector& xlow, const SimpleVector& ixlow, const SimpleVector& xupp, const SimpleVector& ixupp) const;
   public:

//      bool combineColAdaptParent();

      bool reductionsEmpty();

public:
      // todo getter, setter for element access of nnz counter???
      double getObjOffset() const { return objOffset; };
      double addObjOffset(double addOffset);
      void setObjOffset(double offset);
      int getNChildren() const { return nChildren; };

      /// synchronizing the problem over all mpi processes if necessary
      void allreduceLinkingVarBounds();
      void allreduceAndApplyLinkingRowActivities();
      void allreduceAndApplyNnzChanges();
      void allreduceAndApplyBoundChanges();

      /// interface methods called from the presolvers when they detect a possible modification
// todo make bool and ginve feedback or even better - return some enumn maybe?
      void fixColumn(int node, int col, double value);
      bool rowPropagatedBounds( SystemType system_type, int node, BlockType block_type, int row, int col, double ubx, double lbx);
      void removeRedundantRow(SystemType system_type, int node, int row, bool linking);
      void removeParallelRow(SystemType system_type, int node, int row, bool linking);

      /* call whenever a single entry has been deleted from the matrix */
      void deleteEntry(SystemType system_type, int node, BlockType block_type, int row_index, int& index_k, int& row_end);
      void adjustMatrixBoundsBy(SystemType system_type, int node, BlockType block_type, int row_index, double value);
      void updateTransposedSubmatrix( SystemType system_type, int node, BlockType block_type, std::vector<std::pair<int, int> >& elements);

      /* methods for verifying state of presData or querying the problem */
public :
      bool verifyNnzcounters();
      bool elementsDeletedInTransposed() { return elements_deleted == elements_deleted_transposed; };

      bool nodeIsDummy(int node, SystemType system_type) const;
      bool hasLinking(SystemType system_type) const;

private:
/// methods for modifying the problem
      void adjustRowActivityFromDeletion(SystemType system_type, int node, BlockType block_type, int row, int col, double coeff);
      void setVarBoundTo(SystemType system_type, int node, double value); // todo use and implement
      void removeColumn(int node, int col, double fixation);
      void removeColumnFromMatrix(SystemType system_type, int node, BlockType block_type, int col, double fixation);
      void removeRow(SystemType system_type, int node, int row, bool linking);
      void removeRowFromMatrix(SystemType system_type, int node, BlockType block_type, int row);
      void removeEntryInDynamicStorage(SparseStorageDynamic& storage, int row_idx, int col_idx) const;

      void removeIndexRow(SystemType system_type, int node, BlockType block_type, int row_index, int amount);
      void removeIndexColumn(int node, BlockType block_type, int col_index, int amount);
/// methods for querying the problem in order to get certain structures etc.
      SparseGenMatrix* getSparseGenMatrix(SystemType system_type, int node, BlockType block_type);
};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_ */

