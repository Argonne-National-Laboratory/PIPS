/*
 * Scaler.h
 *
 *  Created on: 19.12.2017
 *      Authors: Daniel Rehfeldt, Svenja Uslu
 */

#ifndef PIPS_IPM_CORE_ABSTRACT_SCALER_H_
#define PIPS_IPM_CORE_ABSTRACT_SCALER_H_


class Data;
class Variables;
class Residuals;
class LinearSystem;
class Status;
class OoqpMonitor;
class OoqpStartStrategy;
class ProblemFormulation;

/**  * @defgroup Scaler
 *
 * Interior-point scalers
 * @{
 */

/**
 * Abstract base class for scalers.
 */

enum ScalerType {SCALER_NONE, SCALER_EQUI, SCALER_GEO};

class Scaler
{

public:

  Scaler(Data * prob);
  virtual ~Scaler();

  /** scale */
  virtual void scale( ProblemFormulation * formulation,
            Data * prob, Variables * vars, Residuals * resid) = 0;

  /** unscale */
  virtual void unscale( ProblemFormulation * formulation,
            Data * prob, Variables * vars, Residuals * resid) = 0;

};

//@}


#endif /* PIPS_IPM_CORE_ABSTRACT_SCALER_H_ */
