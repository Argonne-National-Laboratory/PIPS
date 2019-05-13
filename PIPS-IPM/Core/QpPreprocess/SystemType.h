/*
 * SystemType.h
 *
 *  Created on: 02.05.2019
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_SYSTEMTYPE_H_
#define PIPS_IPM_CORE_QPPREPROCESS_SYSTEMTYPE_H_

enum SystemType
{
   EQUALITY_SYSTEM,
   INEQUALITY_SYSTEM
};

enum BlockType
{
   LINKING_VARS_BLOCK,
   CHILD_BLOCK,
   LINKING_CONS_BLOCK
};

#endif /* SYSTEMTYPE_H_ */
