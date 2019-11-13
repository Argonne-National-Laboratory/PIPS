/*
 * StochPostsolver.h
 *
 *  Created on: 02.05.2019
 *      Author: Nils-Christian Kempke
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPOSTSOLVER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPOSTSOLVER_H_

#include <vector>

#include "QpPostsolver.h"
#include "StochVector.h"
#include "sData.h"
#include "sVars.h"
#include "SystemType.h"


class StochPostsolver : public QpPostsolver {
public:

      StochPostsolver( const sData& original_problem );
      virtual ~StochPostsolver();


      void notifySingletonEqualityRow( int node, int row, BlockType block_type, int col, double coeff, double rhs);
      void notifySingletonIneqalityRow( int node, int row, BlockType block_type, int col, double coeff, double lhs, double rhs );


      void notifyRedundantRow( SystemType system_type, int node, unsigned int row, bool linking_constraint, const std::vector<int>& indices_row,
         const std::vector<double> values_row );
      void notifyFixedColumn( int node, unsigned int col, double value, const std::vector<int>& indices_col, const std::vector<double>& values_col);
      void notifyFixedEmptyColumn( int node, unsigned int col, double value);

      void notifyRowPropagated( SystemType system_type, int node, int row, bool linking_constraint, int column, double lb, double ub, double* values, int* indices, int length);
      void notifyDeletedRow( SystemType system_type, int node, int row, bool linking_constraint);
      void notifyParallelColumns();
      void notifyParallelRowSubstitution(SystemType system_type, int node_row, int var1, int row1, int node_var1, int var2, int row2, 
         int node_var2, double scalar, double translation);


      virtual PostsolveStatus postsolve(const Variables& reduced_solution, Variables& original_solution) const;
private:

      int my_rank;
      bool distributed;

      /* can represent a column or row of the problem - EQUALITY/INEQUALITY system has to be stored somewhere else */
      struct INDEX
      {
         INDEX(int node, int index) : node(node), index(index) {};
         int node;
         int index;
      } ;

      enum ReductionType
      {
         FIXED_COLUMN = 0,
         SUBSTITUTED_COLUMN = 1,
         PARALLEL_COLUMN = 2,
         DELETED_ROW = 3,
         REDUNDANT_ROW = 4,
         BOUNDS_TIGHTENED = 5,
         SINGLETON_EQUALITY_ROW = 6,
         SINGLETON_INEQUALITY_ROW = 7,
         FIXED_EMPTY_COLUMN = 8,
         PARALLEL_ROW_SUBSTITUTION = 9,
      };

      const unsigned int n_rows_original;
      const unsigned int n_cols_original;

      /// for now mapping will contain a dummy value for columns that have not been fixed and the value the columns has been fixed to otherwise
      StochVectorBase<int>* padding_origcol;
      StochVectorBase<int>* padding_origrow_equality;
      StochVectorBase<int>* padding_origrow_inequality;

      std::vector<ReductionType> reductions;
      std::vector<INDEX> indices;
      
      std::vector<double> values;
      std::vector<unsigned int> start_idx_values;


      // todo KKTchecker

      void finishNotify();

/// postsolve operations
      void setOriginalVarsFromReduced(const sVars& reduced_vars, sVars& original_vars) const;

      template <typename T>
      void setOriginalValuesFromReduced(StochVectorBase<T>& original_vector,
         const StochVectorBase<T>& reduced_vector,
         const StochVectorBase<int>& padding_original) const;

      template <typename T>
      void setOriginalValuesFromReduced(SimpleVectorBase<T>& original_vector,
         const SimpleVectorBase<T>& reduced_vector,
         const SimpleVectorBase<int>& padding_original) const;


};





#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPOSTSOLVER_H_ */
