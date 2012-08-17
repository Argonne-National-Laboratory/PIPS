/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef MA86QPGENFACTORY
#define MA86QPGENFACTORY

#include "QpGenSparseSeq.h"

class QpGenSparseMa86 : public QpGenSparseSeq {
public:
  QpGenSparseMa86( int nx, int my, int mz,
		       int nnzQ, int nnzA, int nnzC );
  LinearSystem * makeLinsys( Data * prob_in );
};

#endif
