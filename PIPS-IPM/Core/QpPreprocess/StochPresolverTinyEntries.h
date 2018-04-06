/*
 * StochPresolverTinyEntries.h
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERTINYENTRIES_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERTINYENTRIES_H_


class StochPresolverTinyEntries : public StochPresolverBase
{
      StochPresolverTinyEntries(PresolveData& presData);

      ~StochPresolverTinyEntries();

      // remove small matrix entries and return number of eliminations
      virtual bool applyPresolving(int& nelims);

   private:

      int removeTinyEntriesSystemA();
      int removeTinyEntriesSystemC();

      int removeTinyChild( int it, SystemType system_type );

      int removeTinyInnerLoop( int it, SystemType system_type, BlockType block_type );
      void storeRemovedEntryIndex(int rowidx, int colidx, int it, BlockType block_type);

      // data
};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERTINYENTRIES_H_ */
