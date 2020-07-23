/*
 * StochPresolver.h
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_

#include "QpPresolver.h"
#include <vector>

class Data;
class Postsolver;
class StochPresolverBase;
class PresolveData;

/**  * @defgroup QpPreprocess
 *
 * QP presolver
 * @{
 */

/**
 * Derived class for QP presolvers.
 */
class StochPresolver : public QpPresolver
{
private:
   const int my_rank;
   /** limit for max rounds to apply all presolvers */
   const int limit_max_rounds;
   /** should free variables' bounds be reset after presolve (given the row implying these bounds was not removed */
   const bool reset_free_variables_after_presolve;
   /** should the problem be written to std::cout before and after presolve */
   const bool print_problem;
   /** should the presolved problem be written out in MPS format */
   const bool write_presolved_problem;

   const int verbosity;

   PresolveData* presData;
   std::vector<StochPresolverBase*> presolvers;

   void resetFreeVariables();
public:

   StochPresolver(const Data* prob, Postsolver* postsolver);
   virtual ~StochPresolver();

   Data* presolve() override;
};

//@}




#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_ */
