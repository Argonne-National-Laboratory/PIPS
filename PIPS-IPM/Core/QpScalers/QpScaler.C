/*
 * QpScaler.C
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

#include "QpScaler.h"

QpScaler::QpScaler()
{
   vec_rowscaleA = NULL;
   vec_colscaleA = NULL;
   vec_rowscaleC = NULL;
   vec_colscaleC = NULL;
}

QpScaler::~QpScaler()
{
   delete vec_rowscaleA;
   delete vec_colscaleA;
   delete vec_rowscaleC;
   delete vec_colscaleC;
}

void QpScaler::scale( ProblemFormulation * formulation,
            QpGenData * prob, QpGenVars * vars, QpGenResiduals * resid)
{

}

void QpScaler::unscale( ProblemFormulation * formulation,
      QpGenData * prob, QpGenVars * vars, QpGenResiduals * resid)
{

}

