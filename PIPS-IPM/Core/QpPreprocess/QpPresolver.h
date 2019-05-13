/*
 * QpPresolver.h
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_QPPRESOLVER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_QPPRESOLVER_H_

#include "Presolver.h"

class Data;
class Postsolver;
/**  * @defgroup QpPreprocess
 *
 * QP scaler
 * @{
 */

/**
 * Abstract base class for QP Presolvers.
 */

class QpPresolver : public Presolver
{
   protected:
      const Data* const origprob;
      Postsolver* const postsolver;

   public:
      QpPresolver(const Data* prob, Postsolver* postsolver = NULL);
      virtual ~QpPresolver();

      virtual Data* presolve() = 0;
};


//@}


#endif /* PIPS_IPM_CORE_QPPREPROCESS_QPPRESOLVER_H_ */
