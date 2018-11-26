/*
 * Scaler.h
 *
 *  Created on: 19.12.2017
 *      Authors: Daniel Rehfeldt, Svenja Uslu
 */

#ifndef PIPS_IPM_CORE_ABSTRACT_SCALER_H_
#define PIPS_IPM_CORE_ABSTRACT_SCALER_H_
#include "OoqpVector.h"

class Data;

/**  * @defgroup Preprocessing
 *
 * Interior-point scalers
 * @{
 */

/**
 * Abstract base class for scalers.
 */


class Scaler
{
protected:
  Data* const problem;
  const bool do_bitshifting; // only scale by power of two factors?
  const bool with_sides; // consider lhs/rhs?

public:
  Scaler(Data * prob, bool bitshifting = true, bool usesides = true);
  virtual ~Scaler();

  /** scale */
  virtual void scale() = 0;

  /** unscale given objective value */
  virtual double getOrigObj(double objval) = 0;

  /** compute objective value from given primal solution vector */
  virtual OoqpVector* getOrigObj(const OoqpVector& solprimal) = 0;


};

//@}


#endif /* PIPS_IPM_CORE_ABSTRACT_SCALER_H_ */
