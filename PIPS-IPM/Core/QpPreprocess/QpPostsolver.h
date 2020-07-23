/*
 * QpPostsolver.h
 *
 *  Created on: 03.05.2019
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_QPPOSTSOLVER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_QPPOSTSOLVER_H_

#include "Postsolver.h"

class Data;

/**
 * Abstract base class for QP Postsolvers.
 */

class QpPostsolver : public Postsolver
{
   protected:
      const Data& original_problem;

   public:
      QpPostsolver(const Data& prob);
      virtual ~QpPostsolver();

      /** postsolve reduced solution and set original solution accordingly */
      virtual PostsolveStatus postsolve(const Variables& reduced_solution, Variables& original_solution) = 0;
};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_QPPOSTSOLVER_H_ */
