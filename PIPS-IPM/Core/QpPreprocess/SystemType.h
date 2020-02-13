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
   A_MAT,
   B_MAT,
   BL_MAT
};

/* can point to a column or row of the problem - EQUALITY/INEQUALITY system has to be stored somewhere else */
enum IndexType {COL, ROW, EMPTY_INDEX};

struct INDEX
{
   INDEX() : index_type(EMPTY_INDEX), node(-2), index(-1), linking(false), system_type(EQUALITY_SYSTEM){};
   INDEX(IndexType index_type, int node, int index, bool linking = false, SystemType system_type = EQUALITY_SYSTEM) :
      index_type(index_type), node(node), index(index), linking(linking), system_type(system_type){};

   inline bool isRow() const { return index_type == ROW; };
   inline bool isCol() const { return index_type == COL; };

   const IndexType index_type;
   const int node;
   const int index;
   const bool linking;
   const SystemType system_type;
} ;

#endif /* SYSTEMTYPE_H_ */
