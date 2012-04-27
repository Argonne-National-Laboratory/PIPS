/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/memarray.cxx

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



#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "memarray.hxx"

 
using namespace CH_Tools;

namespace CH_Matrix_Classes {

  std::ostream* materrout=&std::cout;        //global error channel

Memarray* Memarrayuser::memarray;

// **************************************************************************
//                                mat_randgen
// **************************************************************************

GB_rand mat_randgen;

// **************************************************************************
//                                MEmessage
// **************************************************************************

int MEmessage(const MatrixError& me)
{
  if (materrout) {
    (*materrout)<<"MatrixError(";
    switch (me.code){
    case ME_unspec: (*materrout)<<"unspecified"; break;
    case ME_range : (*materrout)<<"range"; break;
    case ME_mem   : (*materrout)<<"memory"; break;
    case ME_dim   : (*materrout)<<"dimension"; break;
    case ME_num   : (*materrout)<<"numeric"; break;
    case ME_warning: (*materrout)<<"warning"; break;
    default : (*materrout)<<"?"; break;
    }
    (*materrout)<<",";
    switch (me.mtype){
    case MTindexmatrix: (*materrout)<<"Indexmatrix";break;
    case MTmatrix     : (*materrout)<<"Matrix";break;       
    case MTsymmetric  : (*materrout)<<"Symmatrix";break;
    case MTsparse     : (*materrout)<<"Sparsemat";break;
    case MTsparsesym  : (*materrout)<<"Sparsesym";break;
    default : (*materrout)<<"?"; break;
    }
    (*materrout)<<"):";
    (*materrout)<<me.message<<std::endl;
  }
  if (me.code!=ME_warning) {
    //exit(1);
    abort();
  }
  return 1;
}

//*************************************************************************
//                              size_index
//*************************************************************************

int Memarray::size_index(register long size)
{
 register long i=max_sizes-1;
 size>>=3;
 while((size>>=1)&&(--i)){}
 return int(max_sizes-1-i);
}

//*************************************************************************
//                              index_size
//*************************************************************************

long Memarray::index_size(register int index)
{
 return 32<<index;
}

//*************************************************************************
//                              addr_index
//*************************************************************************

int Memarray::addr_index(register const char* addr)
{
 register unsigned long val=(unsigned long) addr;
 val>>=6;
 return int(val & addr_mask);
}

//*************************************************************************
//                          get_more_entries
//*************************************************************************

int Memarray::get_more_entries()
{
#if (CONICBUNDLE_DEBUG>=70)
 if (materrout) (*materrout)<<"MORE ENTRIES:"<<2*max_entries<<std::endl;
#endif
 Entry *new_store=new Entry[2*max_entries];

 //---- first copy old contents to new entries and let the 
 //     the old entries point to their new replacements 

 register Entry *pold=entry_store;
 register Entry *pnew=new_store;
 register long i;
 for(i=0;i<max_entries;i++,pold++,pnew++){
     pnew->next=pold->next; pold->next=pnew;
     pnew->addr=pold->addr; pold->addr=0;
     pnew->size=pold->size;
     pnew->index=pold->index;
 }

 //---- now redirect next pointers to the new entries
 
 pold=entry_store;
 pnew=new_store;
 for(i=0;i<max_entries;i++,pold++,pnew++){
     if (pnew->next) pnew->next=pnew->next->next; 
 }
 for(i=0;i<max_sizes;i++){
     if (first_free[i].next) first_free[i].next=first_free[i].next->next;
 }
 for(i=0;i<max_addr_entr;i++){
     if (first_used[i].next) first_used[i].next=first_used[i].next->next;
 }

 //---- get rid of old storage and initialize new empty elements

 delete[] entry_store;
 entry_store=new_store;
 pnew=entry_store+max_entries;
 i=max_entries;
 first_empty.next=pnew;
 max_entries*=2;
 for(;i<max_entries-1;i++,pnew++) pnew->next=pnew+1;
 pnew->next=0;

 return 0;
}

//*************************************************************************
//                           Constructor
//*************************************************************************

Memarray::Memarray(long nre,int nrs,int nrbits)
{
  memarray_users=0;
 in_use=0;
 max_entries=nre;
 max_sizes=nrs;
 max_addr_entr=(1<<(unsigned long)nrbits);
 addr_mask=(unsigned long)max_addr_entr-1;
 entry_store=new Entry[max_entries];
 first_free=new Entry[max_sizes];
 first_used=new Entry[max_addr_entr];
 if ((entry_store==0)||(first_free==0)||(first_used==0)){
   MEmessage(MatrixError(ME_unspec,"Memarray-constructor: memory allocation failed",MTglobalfun));
 }
 Entry* ep=entry_store;
 first_empty.next=ep;
 int i=0;
 for(;i<max_entries-1;i++,ep++) ep->next=ep+1;
 ep->next=0;
 for(i=0;i<max_sizes;i++) {
     first_free[i].next=0;
 }
 for(i=0;i<max_addr_entr;i++){
     first_used[i].next=0;
 }
}

//*************************************************************************
//                             ~Memarray 
//*************************************************************************

Memarray::~Memarray()
{
#if (CONICBUNDLE_DEBUG>=1)
  int cnt=0;
  for(int i=0;i<max_addr_entr;i++){
    Entry *used=first_used[i].next;
    while(used) {
      cnt++;
      used=used->next;
    }
  }
  if (cnt>0){
    if (materrout) (*materrout)<<"**** ERROR: Memarray::~Memarray(): destructing memarray while "<<cnt<<" blocks still in use"<<std::endl;
  }
#endif
  delete[] entry_store; entry_store=0;
  delete[] first_free; first_free=0;
  delete[] first_used; first_used=0;
}

//*************************************************************************
//                             get
//*************************************************************************

long Memarray::get(register long size,char*& addr)
{
 addr=0;
 if (size<=0) return 0;
 register Entry* ep=0;
 int si=size_index(size);

 //---- try to find smallest feasible free entry within this class

 register Entry* ip=first_free+si;
 register Entry* ipn;
 while((ipn=ip->next)){
     if (ipn->size>=size){
         ep=ipn;
         ip->next=ep->next;
#if (CONICBUNDLE_DEBUG>=70)
         if (materrout) (*materrout)<<"DA++  "<<setw(5)<<ep->size<<" ("<<setw(5)<<size<<")"<<std::endl;
#endif
         break;
     }
     ip=ipn;
 }

 //---- if there is none, create a new one

 if (ep==0){   //if there was no free entry of this size
     if (first_empty.next==0) {      //if all entries are taken already
         get_more_entries();
     }
     ep=first_empty.next;
     first_empty.next=ep->next;
     long roundupsize=index_size(si);
     ep->size=(roundupsize>size)?roundupsize:size; 
     ep->addr=new char[ep->size];
     if (ep->addr==0) ep->size=0; //allocation not successful
     ep->index=si;
     in_use+=ep->size;
#if (CONICBUNDLE_DEBUG>=70)
     if (materrout) (*materrout)<<"DA==  "<<setw(5)<<ep->size<<std::endl;
#endif
 }  

 //---- put the entry in the front of the taken list  

 si=addr_index(ep->addr);
 ep->next=first_used[si].next;
 first_used[si].next=ep;

 //---- output address and actual size
 
 addr=ep->addr;
 return ep->size;
}

//*************************************************************************
//                             freedarray
//*************************************************************************

int Memarray::free(register void *addr)
{
 if (addr==0) return 0;

 //---- scan list of taken entries for this address

 int ai=addr_index((const char*) addr);
 register Entry *ip=first_used+ai;
 register Entry *ipn;
 while((ipn=ip->next)){
     if (ipn->addr==addr) break;
     ip=ipn;
 }
 if (!ipn) return 1;  //the address is not in the list 
 
 //---- take the entry from the list

 register Entry *ep=0;
 ep=ipn;
 ip->next=ep->next;
 
 //---- compute size-index and link into sorted position

 register long size=ep->size;
 ip=first_free+ep->index;
 while((ipn=ip->next)&&(ipn->size<size)) ip=ipn;
 ep->next=ipn;
 ip->next=ep;
#if (CONICBUNDLE_DEBUG>=70)
 if (materrout) (*materrout)<<"DA--  "<<setw(5)<<ep->size<<std::endl;
#endif
 return 0;
}

}

