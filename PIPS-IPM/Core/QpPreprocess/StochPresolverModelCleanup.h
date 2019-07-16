/*
 * StochPresolverModelCleanup.h
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERMODELCLEANUP_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERMODELCLEANUP_H_

#include "StochPresolverBase.h"

class StochPresolverModelCleanup : public StochPresolverBase
{
   public:
      StochPresolverModelCleanup(PresolveData& presData, const sData& origProb);

      virtual ~StochPresolverModelCleanup();

      // remove small matrix entries
      virtual void applyPresolving();

   private:
      int removed_entries_total;
      int removed_rows_total;

      int removeRedundantRows(SystemType system_type);
      int removeRedundantRows(SystemType system_type, int node);
      int removeRedundantRows(SystemType system_type, int node, bool linking);
      int removeTinyEntriesFromSystem(SystemType system_type);
      int removeTinyInnerLoop(SystemType system_type, int node, BlockType block_type);
      void fixEmptyColumns();
};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERMODELCLEANUP_H_ */
