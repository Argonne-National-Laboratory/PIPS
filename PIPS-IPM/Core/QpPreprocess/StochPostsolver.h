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


class StochPostsolver : public QpPostsolver{
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

      virtual PostsolveStatus postsolve(const Variables& reduced_solution, Variables& original_solution) const;
protected:

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
      };

      const unsigned int n_rows_original;
      const unsigned int n_cols_original;

      /// for now mapping will contain a dummy value for columns that have not been fixed and the value the columns has been fixed to otherwise
      StochVector* padding_origcol;
      StochVector* padding_origrow_equality;
      StochVector* padding_origrow_inequality;

      std::vector<ReductionType> reductions;
      std::vector<INDEX> indices;
      
      std::vector<double> values;
      std::vector<unsigned int> start_idx_values;


      // todo KKTchecker

//      void setReducedValuesInOrigVector(const StochVector& reduced_vector, StochVector& original_vector, const StochVector& mapping_to_original) const;
//      void setReducedValuesInOrigVector(const SimpleVector& reduced_vector, SimpleVector& original_vector, const SimpleVector& mapping_to_original) const;

private:
      void finishNotify();

      SimpleVector& getSimpleVecRowFromStochVec(const OoqpVector& ooqpvec, int node, BlockType block_type) const
         { return getSimpleVecRowFromStochVec(dynamic_cast<const StochVector&>(ooqpvec), node, block_type); };
      SimpleVector& getSimpleVecColFromStochVec(const OoqpVector& ooqpvec, int node) const
         { return getSimpleVecColFromStochVec(dynamic_cast<const StochVector&>(ooqpvec), node); };
      SimpleVector& getSimpleVecRowFromStochVec(const StochVector& stochvec, int node, BlockType block_type) const;
      SimpleVector& getSimpleVecColFromStochVec(const StochVector& stochvec, int node) const;

/// postsolve operations

      void setOriginalVarsFromReduced(const sVars& reduced_vars, sVars& original_vars) const;

      void setOriginalValuesFromReduced(StochVector& original_vector, const StochVector& reduced_vector, const StochVector& padding_original) const;
      void setOriginalValuesFromReduced(SimpleVector& original_vector, const SimpleVector& reduced_vector, const SimpleVector& padding_original) const;


};





#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPOSTSOLVER_H_ */
