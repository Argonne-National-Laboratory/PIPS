/*
 * StochVectorUtilities.h
 *
 *  Created on: 13.11.2019
 *      Author: bzfkempk
 */
#ifndef PIPS_IPM_CORE_STOCHLINEARALGEBRA_STOCHVECTORUTILITIES_H_
#define PIPS_IPM_CORE_STOCHLINEARALGEBRA_STOCHVECTORUTILITIES_H_

#include "StochVector.h"
#include "SimpleVector.h"

#include <string>
#include <iostream>
#include <limits>


/* Utility functions to return for a given StochVectorBase<T> a certain SimpleVectorBase<T> 
 *	either from one of the children (node = 0, ..., nChildren - 1 )
 * or from the parent node = -1
 *
 * for the col vec function SimpleVectorBase<T> vec is returned 
 *
 * for row vec function either the associated SimpleVectorBase<T> vec (linking = false) or vecl (linking = true) is returned
 *
 * asserts existence of these vectors as well as the specified child
 */
template <typename T>
inline SimpleVectorBase<T>& getSimpleVecFromStochVec(const StochVectorBase<T>& stochvec, int node, bool linking)
{
   assert(-1 <= node && node < static_cast<int>(stochvec.children.size()) );

   if(node == -1)
   {
      if( linking )
      {
         assert(stochvec.vecl);
         return dynamic_cast<SimpleVectorBase<T>&>(*(stochvec.vecl));
      }
      else
      {
         assert(stochvec.vec);
         return dynamic_cast<SimpleVectorBase<T>&>(*(stochvec.vec));
      }
   }
   else
   {
      if( linking )
      {
         assert(stochvec.vecl);
         return dynamic_cast<SimpleVectorBase<T>&>(*(stochvec.vecl));
      }
      else
      {
         assert(stochvec.children[node]->vec);
         return dynamic_cast<SimpleVectorBase<T>&>(*(stochvec.children[node]->vec));
      }
   }
}

template <typename T>
inline SimpleVectorBase<T>& getSimpleVecFromRowStochVec(const OoqpVectorBase<T>& ooqpvec, int node, bool linking)
{ return getSimpleVecFromStochVec(dynamic_cast<const StochVectorBase<T>&>(ooqpvec), node, linking); }

template <typename T>
inline SimpleVectorBase<T>& getSimpleVecFromColStochVec(const OoqpVectorBase<T>& ooqpvec, int node)
   { return getSimpleVecFromStochVec(dynamic_cast<const StochVectorBase<T>&>(ooqpvec), node, false); }

template <typename T>
inline T& getSimpleVecFromRowStochVec(const SmartPointer<OoqpVectorBase<T> >& ooqpvec_handle, const INDEX& row)
{
   assert(row.isRow());
   const OoqpVectorBase<T>& ooqp_vec = *ooqpvec_handle;
   return getSimpleVecFromRowStochVec(ooqp_vec, row);
}

template <typename T>
inline T& getSimpleVecFromRowStochVec(const OoqpVectorBase<T>& ooqpvec, const INDEX& row)
{
   assert(row.isRow());
   SimpleVectorBase<T>& vec = getSimpleVecFromStochVec(dynamic_cast<const StochVectorBase<T>&>(ooqpvec), row.getNode(), row.getLinking());
   const int index = row.getIndex();
   assert(0 <= index && index < vec.n);

   return vec[index];
}

template <typename T>
inline T& getSimpleVecFromColStochVec(const SmartPointer<OoqpVectorBase<T> >& ooqpvec_handle, const INDEX& col)
{
   assert(col.isCol());
   const OoqpVectorBase<T>& ooqp_vec = *ooqpvec_handle;
   return getSimpleVecFromColStochVec(ooqp_vec, col);
}

template <typename T>
inline T& getSimpleVecFromColStochVec(const OoqpVectorBase<T>& ooqpvec, const INDEX& col)
{
   assert(col.isCol());
   SimpleVectorBase<T>& vec = getSimpleVecFromStochVec(dynamic_cast<const StochVectorBase<T>&>(ooqpvec), col.getNode(), false);
   const int index = col.getIndex();
   assert(0 <= index && index < vec.n);

   return vec[index];
}

/// clone the structure of a StochVectorBase<T> into one of type U
template<typename T, typename U>
inline StochVectorBase<U>* cloneStochVector(const StochVectorBase<T>& svec)
{
  StochVectorBase<U>* clone;
  
  if( svec.isKindOf(kStochDummy) )
    return new StochDummyVectorBase<U>();

  if( svec.vecl )
    clone = new StochVectorBase<U>( svec.vec->length(), svec.vecl->length(), svec.mpiComm, -1);
  else
    clone = new StochVectorBase<U>( svec.vec->length(), svec.mpiComm);

  for(size_t it = 0; it < svec.children.size(); it++) {
    clone->AddChild( cloneStochVector<T,U>(*svec.children[it]) );
  }
  return clone;
}

template<typename T, typename U>
inline StochDummyVectorBase<U>* cloneStochVector(const StochDummyVectorBase<T>& dummyvec)
{
  return new StochDummyVectorBase<U>();
}

template<typename T, typename U>
inline StochVectorBase<U>* cloneStochVector(const OoqpVectorBase<T>& ooqpvec)
{
  if( ooqpvec.isKindOf(kStochDummy) )
    return cloneStochVector<T, U>( dynamic_cast<const StochDummyVectorBase<T>&>(ooqpvec) );
  else
    return cloneStochVector<T, U>( dynamic_cast<const StochVectorBase<T>&>(ooqpvec) );
}

inline void writeLBltXltUBToStreamAllStringStream(std::stringstream& sout, 
   const SimpleVector& lb, const SimpleVector& ixlow,
   const SimpleVector& ub, const SimpleVector& ixupp)
{
   assert(lb.length() == ixlow.length());
   assert(lb.length() == ub.length());
   assert(lb.length() == ixupp.length());
   
   for( int i = 0; i < lb.length(); i++ )
   {
      const double lbx = PIPSisZero(ixlow[i]) ? -std::numeric_limits<double>::infinity() : lb[i];
      const double ubx = PIPSisZero(ixupp[i]) ? std::numeric_limits<double>::infinity() : ub[i];
   
      sout << lbx << "\t<= x <=\t" << ubx << "\n";
   }
}

inline void writeLBltXltUBtoStringStreamDenseChild( std::stringstream& sout, 
   const StochVector& lb, const StochVector& ixlow,
   const StochVector& ub, const StochVector& ixupp )
{
   assert(lb.children.size() == ixlow.children.size());
   assert(lb.children.size() == ub.children.size());
   assert(lb.children.size() == ixupp.children.size());
   
   if(lb.isKindOf(kStochDummy))
   {
      assert(ixlow.isKindOf(kStochDummy));
      assert(ub.isKindOf(kStochDummy)); assert(ixupp.isKindOf(kStochDummy));
      return;
   }
   
   sout << "--" << std::endl;
   assert(lb.vec); assert(ixlow.vec); assert(ub.vec); assert(ixupp.vec);
   writeLBltXltUBToStreamAllStringStream(sout, 
      dynamic_cast<const SimpleVector&>(*lb.vec),
      dynamic_cast<const SimpleVector&>(*ixlow.vec),
      dynamic_cast<const SimpleVector&>(*ub.vec),
      dynamic_cast<const SimpleVector&>(*ixupp.vec));
      
   
   for( size_t it = 0; it < lb.children.size(); it++ ){
      sout << "-- " << std::endl;
      writeLBltXltUBtoStringStreamDenseChild(sout, *lb.children[it],
         *ixlow.children[it], *ub.children[it], *ixupp.children[it]);   }

   if( lb.vecl )
   {
      assert(ixlow.vecl); assert(ub.vecl); assert(ixupp.vecl);
      sout << "---" << std::endl;
      writeLBltXltUBToStreamAllStringStream(sout,
         dynamic_cast<const SimpleVector&>(*lb.vecl),
         dynamic_cast<const SimpleVector&>(*ixlow.vecl),
         dynamic_cast<const SimpleVector&>(*ub.vecl),
         dynamic_cast<const SimpleVector&>(*ixupp.vecl));
   }
}

inline void writeLBltXltUBtoStreamDense( std::ostream& out,
   const StochVector& lb, const StochVector& ixlow,
   const StochVector& ub, const StochVector& ixupp, MPI_Comm mpiComm )
{
   assert(lb.children.size() == ixlow.children.size());
   assert(lb.children.size() == ub.children.size());
   assert(lb.children.size() == ixupp.children.size());
   
   if( lb.isKindOf(kStochDummy) )
   {
      assert(ixlow.isKindOf(kStochDummy));
      assert(ub.isKindOf(kStochDummy)); assert(ixupp.isKindOf(kStochDummy));
      return;
   }
   
   const int my_rank = PIPS_MPIgetRank(mpiComm);
   const int world_size = PIPS_MPIgetSize(mpiComm);
   bool is_distributed = PIPS_MPIgetDistributed(mpiComm);
   
   MPI_Status status;
   int l;
   std::stringstream sout;

   if( my_rank == 0)
   {
      sout << "----" << std::endl;
      assert(lb.vec); assert(ixlow.vec); assert(ub.vec); assert(ixupp.vec);
      writeLBltXltUBToStreamAllStringStream(sout, 
         dynamic_cast<const SimpleVector&>(*lb.vec),
         dynamic_cast<const SimpleVector&>(*ixlow.vec),
         dynamic_cast<const SimpleVector&>(*ub.vec),
         dynamic_cast<const SimpleVector&>(*ixupp.vec));

      for( size_t it = 0; it < lb.children.size(); it++ )
         writeLBltXltUBtoStringStreamDenseChild(sout, *lb.children[it],
            *ixlow.children[it], *ub.children[it], *ixupp.children[it]);
      
      out << sout.str();
      sout.str(std::string());

      for( int p = 1; p < world_size; p++ )
      {
         MPI_Probe(p, p, mpiComm, &status);
         MPI_Get_count(&status, MPI_CHAR, &l);
         char *buf = new char[l];
         MPI_Recv(buf, l, MPI_CHAR, p, p, mpiComm, &status);
         std::string rowPartFromP(buf, l);
         out << rowPartFromP;
         delete[] buf;
      }
      if( lb.vecl )
      {
         assert(ixlow.vecl); assert(ub.vecl); 
         assert(ixupp.vecl);
         sout << "---" << std::endl;
         writeLBltXltUBToStreamAllStringStream(sout, 
            dynamic_cast<const SimpleVector&>(*lb.vecl), 
            dynamic_cast<const SimpleVector&>(*ixlow.vecl),
            dynamic_cast<const SimpleVector&>(*ub.vecl),
            dynamic_cast<const SimpleVector&>(*ixupp.vecl));
      }
      sout << "----" << std::endl;
      out << sout.str();
   }
   else if( is_distributed )
   { // rank != 0
      for( size_t it = 0; it < lb.children.size(); it++ )
         writeLBltXltUBtoStringStreamDenseChild(sout, *lb.children[it],
            *ixlow.children[it], *ub.children[it], *ixupp.children[it]);

      std::string str = sout.str();
      MPI_Ssend(str.c_str(), str.length(), MPI_CHAR, 0, my_rank, mpiComm);
   }

   if( is_distributed == 1 )
      MPI_Barrier(mpiComm);
}

#endif
