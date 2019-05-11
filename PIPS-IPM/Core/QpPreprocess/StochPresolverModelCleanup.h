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
      StochPresolverModelCleanup(PresolveData& presData, const sData& origProb, StochPostsolver* postsolver);

      virtual ~StochPresolverModelCleanup();

      // remove small matrix entries
      virtual void applyPresolving();

   private:

      int removeRedundantRows(SystemType system_type);
      int removeRedundantRowsBlockwise(SystemType system_type, int node);
      int removeRedundantLinkingRows(SystemType system_type);
      void computeLinkingRowActivity(SystemType system_type, double* minActivity, double* maxActivity, int nLinkRows);
      void checkRedundantLinkingRow(SystemType system_type, double* minActivity, double* maxActivity, int nLinkRows, bool* rowIsRedundant);

      int removeTinyEntriesFromSystem(SystemType system_type);

      int removeTinyInnerLoop( int it, SystemType system_type, BlockType block_type );

      void setPointersCurrentMatrix(int it, SystemType system_type, BlockType block_type);


};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERMODELCLEANUP_H_ */
