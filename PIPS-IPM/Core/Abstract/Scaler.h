/*
 * Scaler.h
 *
 *  Created on: 19.12.2017
 *      Authors: Daniel Rehfeldt, Svenja Uslu
 */

#ifndef PIPS_IPM_CORE_ABSTRACT_SCALER_H_
#define PIPS_IPM_CORE_ABSTRACT_SCALER_H_

class Data;

/**  * @defgroup Scaler
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
  const bool do_bitshifting; // only scale by power of two factors?

public:
  Scaler(Data * prob, bool bitshifting = true);
  virtual ~Scaler();

  /** scale */
  virtual void scale(Data * prob) = 0;

  /** unscale */
  virtual void unscale(Data * prob) = 0;

};

//@}


#endif /* PIPS_IPM_CORE_ABSTRACT_SCALER_H_ */
