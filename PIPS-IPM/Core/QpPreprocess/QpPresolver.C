/*
 * QpPresolver.C
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#include <iostream>
#include "QpPresolver.h"


QpPresolver::QpPresolver(const Data* prob)
 : Presolver(prob), origprob(prob)
{
}

QpPresolver::~QpPresolver()
{
}
