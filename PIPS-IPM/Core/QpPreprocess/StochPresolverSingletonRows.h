/*
 * StochPresolverSingletonRows.h
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_

#include "StochPresolverBase.h"
#include "SparseStorageDynamic.h"
#include "PresolveData.h"

class StochPresolverSingletonRows : public StochPresolverBase
{
public:
   StochPresolverSingletonRows(PresolveData& presData, const sData& origProb);

   ~StochPresolverSingletonRows();

   // remove singleton rows
   virtual void applyPresolving();

private:
   long long removed_rows; 

   std::vector<int> buffer_found_singleton_equality;
   std::vector<INDEX> buffer_rows_lower;
   std::vector<INDEX> buffer_rows_upper;

   std::vector<double> buffer_xlows;
   std::vector<double> buffer_xupps;
   std::vector<double> buffer_coeffs_lower;
   std::vector<double> buffer_coeffs_upper;

   bool removeSingletonRow( const INDEX& row );
   void getBoundsAndColFromSingletonRow( const INDEX& row, int& node_col, int& col_idx, double& xlow_new, double& xupp_new, double& coeff);

   void removeSingletonLinkingColsSynced();
   void resetBuffers();
};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_ */
