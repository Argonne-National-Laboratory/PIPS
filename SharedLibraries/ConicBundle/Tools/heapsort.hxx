/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Tools/heapsort.hxx

    This file is part of ConciBundle, a C/C++ library for convex optimization.

    ConicBundle is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ConicBundle is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

***************************************************************************** */



#ifndef CH_TOOLS__HEAPSORT_HXX
#define CH_TOOLS__HEAPSORT_HXX

/**  @file heapsort.hxx
    @brief Header declaring template functions for sorting an array of index-objects indexing an array of value-objects comparable by "<".
    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg
*/


#include "mymath.hxx"

namespace CH_Tools {

  /**@defgroup gheapsort Heapsort
     @brief Template functions for heapsorting an array of index-objects indexing an array of value-objects comparable by "<".

  For given                                                  
    
  -  n :    number of indices to be sorted               
  -  ind :  a class indexable by integers,                  
            a function swap(ind[i],ind[j]) is required    
            which swaps the contents of ind[i] and ind[j] 
  -  val :  class indexable by elements of ind,           
            for two vals val[i] and val[j] operator "<"     
            (val[i]<val[j]) yielding 0/1 must be defined  
                                                             
  the routines sort the first n entries in the index array such that      
                                                             
         val[ind[i]]<=val[ind[j]] for 0<=i<j<n                
  */
  //@{
 
  /**@brief Heapify element i (in 0,...,n-1) in ind, so that afterwards val[ind[j]]>=max(val[ind[2*j+1]],val[ind[2*j+2]] for the subtree with root i
     @param i index of element to be heapfied (indexing starts with 0)
     @param n index array goes from 0 .. n-1 
     @param ind index array into val, rearrangeable by swap(ind[j],ind[k])
     @param val value array index by ind[j], comparable by "<"
  */
  template<class I,class V> void heapify(int i, int n, I& ind,const V& val);

  /**@brief Build a heap out of ind[0] to ind[n-1], so that val[ind[j]]>=max(val[ind[2*j+1]],val[ind[2*j+2]] for all j
     @param n index array goes from 0 .. n-1 
     @param ind index array into val, rearrangeable by swap(ind[j],ind[k])
     @param val value array index by ind[j], comparable by "<"
  */
  template<class I,class V> void build_heap(int n, I& ind, const V& val);

  /**@brief Sort ind[0]...ind[n-1] so that val[ind[i]]<=val[ind[j]] for 0<=i<j<n  
     @param n index array goes from 0 .. n-1 
     @param ind index array into val, rearrangeable by swap(ind[j],ind[k])
     @param val value array index by ind[j], comparable by "<"
  */     
  template<class I,class V> void heapsort(int n, I& ind, const V& val);


  //@}

//--------------------------------------------------------------
//                    heapsort -- implementation
//--------------------------------------------------------------


template<class I,class V> inline void heapify(int i,int n,I& ind,const V& val)
{
 int k,j;
 k=i;
 while(2*k+1<n){
   using namespace CH_Matrix_Classes;
   if (2*k==n-2){
     if (val[ind[k]]<val[ind[2*k+1]])
       swap(ind[k],ind[2*k+1]);
     break;
   }
   if ((val[ind[k]]<val[ind[2*k+1]])||
       (val[ind[k]]<val[ind[2*k+2]])){
     if (val[ind[2*k+1]]>val[ind[2*k+2]]) j=2*k+1;
         else j=2*k+2;
     swap(ind[k],ind[j]);
     k=j;
   }
   else break;
 }
}

template<class I,class V> inline void build_heap(int n,I& ind,const V& val)
{
 for(int i=n/2;--i>=0;) heapify(i,n,ind,val);
}

template<class I,class V> inline void heapsort(int n,I& ind,const V& val)
{
 int i;

 //build heap
 build_heap(n,ind,val);

 //extract greatest and rebuild heap
 for(i=n;--i>=1;){
   using namespace CH_Matrix_Classes;
   swap(ind[i],ind[0]);
   heapify(0,i,ind,val);
 }
}



}

#endif

