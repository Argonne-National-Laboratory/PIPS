/*
 * StochPresolverModelCleanup.h
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERMODELCLEANUP_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERMODELCLEANUP_H_

#include "StochPresolverBase.h"
#include "pipsport.h"

class StochPresolverModelCleanup : public StochPresolverBase
{
   public:
      StochPresolverModelCleanup(PresolveData& presData, const sData& origProb);

      virtual ~StochPresolverModelCleanup();

      // remove small matrix entries
      bool applyPresolving() override;

   private:
      /** limit for the size of a matrix entry below which it will be removed from the problem */
      const double limit_min_mat_entry;
      /** max for the matrix entry when the impact of entry times (bux-blx) is considered */
      const double limit_max_matrix_entry_impact;
      /** difference in orders between feastol and the impact of entry times (bux-blx) for an entry to get removed */
      const double limit_matrix_entry_impact_feasdist;

      int removed_entries_total;
      int fixed_empty_cols_total;
      int removed_rows_total;

      int removeRedundantRows(SystemType system_type);
      int removeRedundantRows(SystemType system_type, int node);
      int removeRedundantRows(SystemType system_type, int node, bool linking);
      int removeTinyEntriesFromSystem(SystemType system_type);
      int removeTinyInnerLoop(SystemType system_type, int node, BlockType block_type);
      int fixEmptyColumns();
};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERMODELCLEANUP_H_ */
