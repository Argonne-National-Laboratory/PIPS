/*
 * Postsolver.h
 *
 *  Created on: 03.05.2019
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_ABSTRACT_POSTSOLVER_H_
#define PIPS_IPM_CORE_ABSTRACT_POSTSOLVER_H_


class Data;
class Variables;

enum PostsolveStatus
{
   PRESOLVE_OK, PRESOLVE_FAIL // TODO
};

/**
 * Abstract base class for Postsolvers.
 */

class Postsolver
{
public:
  Postsolver(const Data& prob) {};
  virtual ~Postsolver() {};

  /** postsolve reduced solution and set original solution accordingly */
  virtual PostsolveStatus postsolve(const Variables& reduced_solution, Variables& original_solution) = 0;

};

//@}


#endif /* PIPS_IPM_CORE_ABSTRACT_POSTSOLVER_H_ */
