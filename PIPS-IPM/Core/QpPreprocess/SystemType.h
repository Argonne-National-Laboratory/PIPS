/*
 * SystemType.h
 *
 *  Created on: 02.05.2019
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_SYSTEMTYPE_H_
#define PIPS_IPM_CORE_QPPREPROCESS_SYSTEMTYPE_H_

#include <ostream>

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
public:
   INDEX() : index_type(EMPTY_INDEX), node(-2), index(-1), linking(false), system_type(EQUALITY_SYSTEM){};

   INDEX(IndexType index_type, int node, int index, bool linking = false, SystemType system_type = EQUALITY_SYSTEM) :
      index_type(index_type), node(node), index(index), linking(linking), system_type(system_type)
   {
      if(linking)
         assert( node == -1 );
      assert( 0 <= index );
      assert( -1 <= node );
   };

   bool operator==(const INDEX& i) const {
      return index_type == i.index_type && node == i.node && index == i.index && linking == i.linking && system_type == i.system_type;
   }

   bool operator!=(const INDEX& i) const {
      return index_type != i.index_type || node != i.node || index != i.index || linking != i.linking || system_type != i.system_type;
   }

   inline bool isRow() const { return index_type == ROW; };
   inline bool isCol() const { return index_type == COL; };
   inline bool isEmpty() const { return index_type == EMPTY_INDEX; };
   inline bool isLinkingCol() const { assert(this->isCol()); return node == -1; };
   inline bool isLinkingRow() const { assert(this->isRow()); return linking; };
   inline bool hasValidNode(int nChildren) const { return -1 <= node && node < nChildren; };
   inline bool inEqSys() const { return system_type == EQUALITY_SYSTEM; }
   inline bool inInEqSys() const { return system_type == INEQUALITY_SYSTEM; }

   BlockType getBlockOfColInRow( const INDEX& col ) const
   {
      assert(isRow());
      assert(col.isCol());

      /* linking row */
      if(linking)
         return BL_MAT;
      else if( col.isLinkingCol() )
      {
         /* row node == -1 and col node == -1 and it is not linking -> B0_Mat */
         if( node == -1 )
            return B_MAT;
         /* row node != -1 and col node == -1 -> some A_MAT */
         else
            return A_MAT;
      }
      /* not a linking row/col and not B0 -> B_MAT */
      else
         return B_MAT;
   }

   friend std::ostream& operator<< (std::ostream &out, const INDEX& i)
   {
      if(i.isRow())
      {
         out << "INDEX(ROW," << i.node << "," << i.index << "," << (i.linking ? "true" : "false") << "," << ((i.system_type == EQUALITY_SYSTEM) ? "EQU_SYS" : "INEQ_SYS") << ")";
      }
      else if(i.isCol())
      {
         out << "INDEX(COL," << i.node << "," << i.index << ")";
      }
      else
      {
         out << "INDEX(EMPTY_INDEX)" << std::endl;
      }
      return out;
   }

   inline const IndexType& getType() const { return index_type; };
   inline const int& getNode() const { return node; };
   inline const int& getIndex() const { return index; };
   inline const bool& getLinking() const { return linking; };
   inline const SystemType& getSystemType() const { return system_type; };

private:
   IndexType index_type;
   int node;
   int index;
   bool linking;
   SystemType system_type;
};


#endif /* SYSTEMTYPE_H_ */
