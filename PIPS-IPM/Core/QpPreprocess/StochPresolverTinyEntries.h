/*
 * StochPresolverTinyEntries.h
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERTINYENTRIES_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERTINYENTRIES_H_

#include "StochPresolverBase.h"
#include <math.h>

class StochPresolverTinyEntries : public StochPresolverBase
{
   public:
      StochPresolverTinyEntries(PresolveData& presData);

      ~StochPresolverTinyEntries();

      // remove small matrix entries
      virtual void applyPresolving();

   private:

      int removeTinyEntriesSystemA();
      int removeTinyEntriesSystemC();

      int removeTinyChild( int it, SystemType system_type );

      int removeTinyInnerLoop( int it, SystemType system_type, BlockType block_type );
      /** initialize current pointer for matrices and vectors.
       * If it==-1, we are at parent and want block B_0 (Bmat).
       * Returns false if it is a dummy child. */
      bool updateCPforTinyEntry(int it, SystemType system_type);

      // data
};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERTINYENTRIES_H_ */
