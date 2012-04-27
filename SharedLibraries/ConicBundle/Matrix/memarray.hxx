/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/memarray.hxx

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



#ifndef CH_MATRIX_CLASSES__MEMARRAY_HXX
#define CH_MATRIX_CLASSES__MEMARRAY_HXX

/**  @file memarray.hxx
    @brief Header for simple memory management tools that support frequent allocation and deallocation of arrays of the same sizes.

    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg

*/ 

#include <iostream>
#include "matop.hxx"
#if (CONICBUNDLE_DEBUG>=1)
#include <iomanip>
#include <stdlib.h>
#endif


namespace CH_Matrix_Classes {

  /**@defgroup memarraygroup Classes and Functions for Memory Management
     @brief implements a simple approach to support frequent allocation and deallocation of arrays of the same sizes.
  */
  //@{

  /** @brief A simple memory manager for frequent allocation and deallocation of arrays of roughly the same size.

      The Manager only allocates blocks of size 2^n and keeps for each n
      a singly linked list holding the free blocks of size 2^n. The occupied
      blocks are kept in singly linked lists indexed by the last few bits of
      their addresses in order to facilitate fast retrieval of the corresponding
      memory management entry. Currently, a block once allocated is not freed
      again unless the memory manager is destructed. 

      Information about the allocated blocks is stored in #Entry items that also
      serve for forming the linked lists. A large array of these items is 
      allocated initially and whenever no more free Entry items are available 
      (i.e. all are filled with allocated blocks) a new array of Entry items,
      double in size, is allocated and the old information of the old array
      is copied to the first half of the new items. 
  */
  class Memarray{
   private:

    ///holds the information of one allocated block and serves as an item in the singly linked lists
    class Entry{
     public:
      Entry *next;   ///< points to the next entry in the list
      char *addr;    ///< points to the allocated block
      long size;     ///< gives the size of the allocated block
      int index;     ///< gives the index of the size group of this block
                   
      ///
      Entry(){next=0; addr=0; size=0;}
      ///
      ~Entry(){delete[] addr;}
    };

    long max_entries;        ///< current number of #Entry items available
    long max_sizes;          ///< current number of lists holding free blocks
    long max_addr_entr;      ///< current number of lists hodling occupied blocks
    unsigned long addr_mask; ///< mask to extract last bits of an address as index for freeing
    unsigned long in_use;    ///< number of #Entry items in use (pointing to an allocated block)
    unsigned long memarray_users; ///< number of objects announced as "living" users of this memory manager 
    
    Entry first_empty;  ///< its next pointer points to the first free #Entry item, that does not yet hold an allocated block 
    Entry* entry_store; ///< points to the allocated array of #Entry items
    Entry* first_free;  ///< points to an array of size max_sizes, its Entry i start the free list for size 2^i
    Entry* first_used;  ///< points to an array of size max_addr_entr, its Entry i start the list of entries whose memory is in use and the last bits of the address form the number i 
    
    /// compute index for free list to a request of size
    int size_index(register long size); 
    /// compute size of blocks in the free list with this index
    long index_size(register int index);
    /// compute index into address class to retrieve entry of "addr"
    int addr_index(register const char *addr); 
    /// double number of available #Entry items
    int get_more_entries();
    
   public:
    /// specify inital number of Entry items, the number of 2^i classes, the number of last bits used in finding entries to used adresses
    Memarray(long number_of_entries,int number_of_sizes,int address_bits);
    ///
    ~Memarray();
    /// returns the number of #Entry items in use
    unsigned long get_in_use() const {return in_use;}
    /// returns the number of announced "living" users of this Memory manager
    unsigned long get_memarray_users() const {return memarray_users;}
    /// announce a new user
    unsigned long increment_memarray_users() {return ++memarray_users;}
    /// decrement number of announced users 
    unsigned long decrement_memarray_users() {return --memarray_users;}
   
    /// request a character array of at least this size, the address is then stored in addr and the actual size is returned  
    long get(register long size,char *& addr);
    /// request a double array of at least this size, the address is then stored in addr and the actual size is returned  
    long get(long size,double *& addr)
    {return get(size*long(sizeof(double)),(char*&)addr)/long(sizeof(double));}
    /// request a float array of at least this size, the address is then stored in addr and the actual size is returned  
    long get(long size,float *& addr)
    {return get(size*long(sizeof(float)),(char*&)addr)/long(sizeof(float));}
    /// request a long array of at least this size, the address is then stored in addr and the actual size is returned  
    long get(long size,long *& addr)
    {return get(size*long(sizeof(long)),(char*&)addr)/long(sizeof(long));}
    /// request an int array of at least this size, the address is then stored in addr and the actual size is returned  
    long get(long size,int *& addr)
    {return get(size*long(sizeof(int)),(char*&)addr)/long(sizeof(int));}
    /// free the array pointed to by addr (addr must be an address returned by get)
    int free(register void *addr);
  };

  /// All derived classes share a common Memarray memory manager, which is generated with the first user and destructed when the last user is destructed. 
  class Memarrayuser
  {
   protected: 
    /// pointer to common memory manager for all Memarrayusers, instantiated in memarray.cxx
    static Memarray* memarray;
   public:
    /// if #memarray is NULL, then a new Memarray is generated. In any case the number of users of the Memarray is incremented
    Memarrayuser()
    { 
      if (memarray==0) {
	memarray=new Memarray(1L,60,10);
      }
      memarray->increment_memarray_users();
    }
    
    ///the number of users is decremented and the Memarray memory manager is destructed, if the number is zero.
    virtual ~Memarrayuser()
    { 
#if (CONICBUNDLE_DEBUG>=1)
      if (memarray==0) {
	MEmessage(MatrixError(ME_unspec,"*** Error: Memarrayuser::~Memarrayuser(): memory management killed prematurely",MTglobalfun));
      }
#endif
      if (memarray->decrement_memarray_users()==0) {
	delete memarray;
	memarray=0;
      }   
    }
  };
  
  /// provide sufficient memory for an existing array, reallocating and copying the old information upon need, returns 0 upon success, !=0 upon failure. 
  template<class T>
  inline int mem_provide(
			 Memarray& memarray,  ///< the Memory Manger which provided the address in store
			 long provide,        ///< the minimum number of entries needed
			 long in_use,         ///< the number of entries currently in use
			 long& avail,         ///< the number of entries currently available and available afterwards
			 T*& store            ///< the current address and the ouput address of the starting element
			 )
  {
    if (provide<avail) return 0;
    T* tmpstore;
    long tmpavail=memarray.get(((provide>2*avail)?provide:2*avail)*long(sizeof(T)),(char *&)tmpstore)/long(sizeof(T));
    if ((tmpstore==0)&&(provide<2*avail)){
      tmpavail=memarray.get(provide*long(sizeof(T)),(char *&)tmpstore)/long(sizeof(T));
    }
    if (tmpstore==0) return 1;
    long i;
    for(i=0;i<in_use;i++){
      tmpstore[i]=store[i];
    }
    memarray.free((void *)store);
    store=tmpstore;
    avail=tmpavail;
    return 0;
  }
  
  /// provide sufficient memory for an existing array, reallocating and copying the old information and initializing the new entries to 0, returns 0 upon success, !=0 upon failure. 
  template<class T>
  inline int mem_provide_init0(
			       Memarray& memarray, ///< the Memory Manger which provided the address in store
			       long provide,       ///< the minimum number of entries needed
			       long& avail,        ///< the number of entries currently available and available afterwards
			       T*& store           ///< the current address and the ouput address of the starting element
			       )
  {
    if (provide<avail) return 0;
    T* tmpstore;
    long tmpavail=memarray.get(((provide>2*avail)?provide:2*avail)*long(sizeof(T)),(char *&)tmpstore)/long(sizeof(T));
    if ((tmpstore==0)&&(provide<2*avail)){
      tmpavail=memarray.get(provide*long(sizeof(T)),(char *&)tmpstore)/long(sizeof(T));
    }
    if (tmpstore==0) return 1;
    long i;
    for(i=0;i<avail;i++){
     tmpstore[i]=store[i];
    }
    for(;i<tmpavail;i++){
      tmpstore[i]=0;
    }
    memarray.free((void *)store);
    store=tmpstore;
    avail=tmpavail;
    return 0;
  }
  
  //@}

}

#endif

