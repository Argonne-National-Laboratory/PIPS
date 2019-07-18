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
class Variables;
class Residuals;

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
  Scaler(Data * prob, bool bitshifting = false, bool usesides = false);
  virtual ~Scaler();

  /** scale */
  virtual void scale() = 0;

  /** unscale given objective value */
  virtual double getOrigObj(double objval) const = 0;

  /** compute original variables from given ones */
  virtual Variables* getUnscaledVariables(const Variables& vars) const = 0;

  /** compute original residuals from given ones */
  virtual Residuals* getUnscaledResiduals(const Residuals& resids) const = 0;

  /** compute original vector from given primal vector */
  virtual OoqpVector* getOrigPrimal(const OoqpVector& solprimal) const = 0;

  /** compute original vector from given dual vector */
  virtual OoqpVector* getOrigDualEq(const OoqpVector& soldual) const = 0;

  /** compute original vector from given dual vector */
  virtual OoqpVector* getOrigDualIneq(const OoqpVector& soldual) const = 0;

  /** compute original vector from given dual vector */
  virtual OoqpVector* getOrigDualVarBoundsUpp(const OoqpVector& soldual) const = 0;

  /** compute original vector from given dual vector */
  virtual OoqpVector* getOrigDualVarBoundsLow(const OoqpVector& soldual) const = 0;
};

//@}


#endif /* PIPS_IPM_CORE_ABSTRACT_SCALER_H_ */
