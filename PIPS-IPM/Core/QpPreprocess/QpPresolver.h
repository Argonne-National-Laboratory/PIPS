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
   public:
      QpPresolver(const Data* prob);
      virtual ~QpPresolver();

      virtual Data* presolve() = 0;
};


//@}


#endif /* PIPS_IPM_CORE_QPPREPROCESS_QPPRESOLVER_H_ */
