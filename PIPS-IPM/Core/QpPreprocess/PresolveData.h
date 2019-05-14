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

typedef struct
{
   int colIdx;
   double val;
} COLUMNFORDELETION;

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

struct col_is_smaller
{
    bool operator()(const COLUMNFORDELETION& x, const COLUMNFORDELETION& y) const
    {
        return x.colIdx < y.colIdx;
    }
};

class PresolveData
{

   public:
      sData* presProb;
private:
      StochPostsolver* const postsolver;

      bool outdated_activities;
      bool outdated_bounds;
      bool outdated_nnzs;
public:
      // todo why handle? ..
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

      StochVectorHandle max_act_eq;
      StochVectorHandle min_act_eq;
      StochVectorHandle max_act_ineq;
      StochVectorHandle min_act_ineq;


      /* stuff for handling the update and changes of activities of certain rows */
      StochVectorHandle actmax_eq;
      StochVectorHandle actmin_eq;

      StochVectorHandle actmax_ineq;
      StochVectorHandle actmin_ineq;

      /* for better MPI communication we allocate one contiguous array in storage and
       * make four SimppleVectors each pointing to a part of it */
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


   private:
      int my_rank;
      bool distributed;

      // number of children
      int nChildren;

      // objective offset created by presolving
      double objOffset;
      double obj_offset_chgs;

      int elements_deleted;
      int elements_deleted_transposed;

      std::vector<COLUMNFORDELETION> linkingVariablesMarkedForDeletion;

      /* methods for initializing the object */
   private:
      // initialize row and column nnz counter
      void initNnzCounter();
      void initSingletons();
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
//      COLUMNFORDELETION getColAdaptParent(int i) const;
//      int getNumberColAdParent() const;
//      void addColToAdaptParent(COLUMNFORDELETION colToAdapt);
//      void clearColAdaptParent();

      void allreduceAndApplyLinkingRowActivities();
      void allreduceAndApplyNnzChanges();
      void allreduceAndApplyBoundChanges();


public:
      /* call whenever a single entry has been deleted from the matrix */
      void deleteEntry(SystemType system_type, int node, BlockType block_type, SparseStorageDynamic* storage,
            int row_index, int& index_k, int& row_end);
      void adjustMatrixBoundsBy(SystemType system_type, int node, BlockType block_type, int row_index, double value);
      void updateTransposedSubmatrix( SparseStorageDynamic* transposed, std::vector<std::pair<int, int> >& elements);

      void removeColumn();
      void removeRedundantRow(SystemType system_type, int node, int row, bool linking);
private:
      void removeRow(SystemType system_type, int node, int row, bool linking);
      void removeRowFromMatrix(SystemType system_type, int node, BlockType block_type, int row);
      void removeEntryInDynamicStorage(SparseStorageDynamic& storage, int row_idx, int col_idx) const;

      /* methods for verifying state of presData */
public :
      bool verifyNnzcounters();
      bool elementsDeletedInTransposed() { return elements_deleted == elements_deleted_transposed; };
private:
      void getSparseGenMatrix(SystemType system_type, int node, BlockType block_type, SparseGenMatrix* mat);
      void removeIndexRow(SystemType system_type, int node, BlockType block_type, int row_index, int amount);
      void removeIndexColumn(int node, BlockType block_type, int col_index, int amount);

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_ */

