/*
 * QpPostsolver.C
 *
 *  Created on: 03.05.2019
 *      Author: bzfkempk
 */

#include "QpPostsolver.h"

QpPostsolver::QpPostsolver(const Data& prob)
 : Postsolver(prob), original_problem(prob)
{
}

QpPostsolver::~QpPostsolver()
{
}
