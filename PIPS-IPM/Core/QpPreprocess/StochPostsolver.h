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

      void notifyFixedColumn( int node, unsigned int col, double value);
      void notifyDeletedRow( SystemType system_type, int node, unsigned int row, bool linking_constraint);
      void notifyParallelColumns();

      virtual PostsolveStatus postsolve(const Variables& reduced_solution, Variables& original_solution) const;
protected:

      /* can represent a column or row of the problem - EQUALITY/INEQUALITY system has to be stored somewhere else */
      typedef struct
      {
         int node;
         int index;
      } INDEX;

      enum ReductionType
      {
         FIXED_COLUMN = 0,
         SUBSTITUTED_COLUMN = 1,
         PARALLEL_COLUMN = 2,
         DELETED_ROW = 3
      };

      const unsigned int n_rows_original;
      const unsigned int n_cols_original;

//      StochVector* mapping_to_origcol;
//      StochVector* mapping_to_origrow_equality;
//      StochVector* mapping_to_origrow_inequality;

      std::vector<ReductionType> reductions;
      std::vector<INDEX> indices;
      std::vector<double> values;
      std::vector<unsigned int> start_idx_values;


      // todo KKTchecker

//      void setReducedValuesInOrigVector(const StochVector& reduced_vector, StochVector& original_vector, const StochVector& mapping_to_original) const;
//      void setReducedValuesInOrigVector(const SimpleVector& reduced_vector, SimpleVector& original_vector, const SimpleVector& mapping_to_original) const;

private:
      void finishNotify();

};





#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPOSTSOLVER_H_ */
