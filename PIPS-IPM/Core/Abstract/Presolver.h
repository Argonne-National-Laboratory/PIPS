/*
 * Presolver.h
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_ABSTRACT_PRESOLVER_H_
#define PIPS_IPM_CORE_ABSTRACT_PRESOLVER_H_


class Data;

/**  * @defgroup Preprocessing
 *
 * Interior-point presolvers
 * @{
 */

/**
 * Abstract base class for presolvers.
 */


class Presolver
{
protected:

public:
  Presolver(const Data * prob);
  virtual ~Presolver();

  /** presolve */
  virtual Data* presolve(const Data* prob) = 0;

};

//@}


#endif /* PIPS_IPM_CORE_ABSTRACT_PRESOLVER_H_ */
