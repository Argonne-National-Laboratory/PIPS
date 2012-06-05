/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   misc.c
 * @brief  miscellaneous methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/misc.h"
#include "scip/intervalarith.h"

#include "scip/struct_misc.h"



/*
 * Priority Queue
 */

#define PQ_PARENT(q) (((q)+1)/2-1)
#define PQ_LEFTCHILD(p) (2*(p)+1)
#define PQ_RIGHTCHILD(p) (2*(p)+2)


/** resizes element memory to hold at least the given number of elements */
static
SCIP_RETCODE pqueueResize(
   SCIP_PQUEUE*          pqueue,             /**< pointer to a priority queue */
   int                   minsize             /**< minimal number of storable elements */
   )
{
   assert(pqueue != NULL);
   
   if( minsize <= pqueue->size )
      return SCIP_OKAY;

   pqueue->size = MAX(minsize, (int)(pqueue->size * pqueue->sizefac));
   SCIP_ALLOC( BMSreallocMemoryArray(&pqueue->slots, pqueue->size) );

   return SCIP_OKAY;
}

/** creates priority queue */
SCIP_RETCODE SCIPpqueueCreate(
   SCIP_PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int                   initsize,           /**< initial number of available element slots */
   SCIP_Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   )
{
   assert(pqueue != NULL);
   assert(ptrcomp != NULL);

   initsize = MAX(1, initsize);
   sizefac = MAX(1.0, sizefac);

   SCIP_ALLOC( BMSallocMemory(pqueue) );
   (*pqueue)->len = 0;
   (*pqueue)->size = 0;
   (*pqueue)->sizefac = sizefac;
   (*pqueue)->slots = NULL;
   (*pqueue)->ptrcomp = ptrcomp;
   SCIP_CALL( pqueueResize(*pqueue, initsize) );

   return SCIP_OKAY;
}

/** frees priority queue, but not the data elements themselves */
void SCIPpqueueFree(
   SCIP_PQUEUE**         pqueue              /**< pointer to a priority queue */
   )
{
   assert(pqueue != NULL);

   BMSfreeMemoryArray(&(*pqueue)->slots);
   BMSfreeMemory(pqueue);
}

/** clears the priority queue, but doesn't free the data elements themselves */
void SCIPpqueueClear(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);

   pqueue->len = 0;
}

/** inserts element into priority queue */
SCIP_RETCODE SCIPpqueueInsert(
   SCIP_PQUEUE*          pqueue,             /**< priority queue */
   void*                 elem                /**< element to be inserted */
   )
{
   int pos;

   assert(pqueue != NULL);
   assert(pqueue->len >= 0);
   assert(elem != NULL);

   SCIP_CALL( pqueueResize(pqueue, pqueue->len+1) );

   /* insert element as leaf in the tree, move it towards the root as long it is better than its parent */
   pos = pqueue->len;
   pqueue->len++;
   while( pos > 0 && (*pqueue->ptrcomp)(elem, pqueue->slots[PQ_PARENT(pos)]) < 0 )
   {
      pqueue->slots[pos] = pqueue->slots[PQ_PARENT(pos)];
      pos = PQ_PARENT(pos);
   }
   pqueue->slots[pos] = elem;

   return SCIP_OKAY;
}

/** removes and returns best element from the priority queue */
void* SCIPpqueueRemove(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   void* root;
   void* last;
   int pos;
   int childpos;
   int brotherpos;

   assert(pqueue != NULL);
   assert(pqueue->len >= 0);
   
   if( pqueue->len == 0 )
      return NULL;

   /* remove root element of the tree, move the better child to its parents position until the last element
    * of the queue could be placed in the empty slot
    */
   root = pqueue->slots[0];
   last = pqueue->slots[pqueue->len-1];
   pqueue->len--;
   pos = 0;
   while( pos <= PQ_PARENT(pqueue->len-1) )
   {
      childpos = PQ_LEFTCHILD(pos);
      brotherpos = PQ_RIGHTCHILD(pos);
      if( brotherpos <= pqueue->len && (*pqueue->ptrcomp)(pqueue->slots[brotherpos], pqueue->slots[childpos]) < 0 )
         childpos = brotherpos;
      if( (*pqueue->ptrcomp)(last, pqueue->slots[childpos]) <= 0 )
         break;
      pqueue->slots[pos] = pqueue->slots[childpos];
      pos = childpos;
   }
   assert(pos <= pqueue->len);
   pqueue->slots[pos] = last;

   return root;
}

/** returns the best element of the queue without removing it */
void* SCIPpqueueFirst(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   if( pqueue->len == 0 )
      return NULL;

   return pqueue->slots[0];
}

/** returns the number of elements in the queue */
int SCIPpqueueNElems(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   return pqueue->len;
}

/** returns the elements of the queue; changing the returned array may destroy the queue's ordering! */
void** SCIPpqueueElems(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   return pqueue->slots;
}




/*
 * Hash Table
 */

/** table of some prime numbers */
static int primetable[] = {
   2,
   7,
   19,
   31,
   59,
   227,
   617,
   1523,
   3547,
   8011,
   17707,
   38723,
   83833,
   180317,
   385897,
   821411,
   1742369,
   3680893,
   5693959,
   7753849,
   9849703,
   11973277,
   14121853,
   17643961,
   24273817,
   32452843,
   49979687,
   67867967,
   86028121,
   104395301,
   122949823,
   141650939,
   160481183,
   179424673,
   198491317,
   217645177,
   256203161,
   314606869,
   373587883,
   433024223,
   492876847,
   553105243,
   613651349,
   694847533,
   756065159,
   817504243,
   879190747,
   941083981,
   982451653,
   INT_MAX
};
static const int primetablesize = sizeof(primetable)/sizeof(int);

/** returns a reasonable hash table size (a prime number) that is at least as large as the specified value */
int SCIPcalcHashtableSize(
   int                   minsize             /**< minimal size of the hash table */
   )
{
   int pos;

   (void) SCIPsortedvecFindInt(primetable, minsize, primetablesize, &pos);
   assert(pos < primetablesize);

   return primetable[pos];
}

/** appends element to the hash list */
static
SCIP_RETCODE hashtablelistAppend(
   SCIP_HASHTABLELIST**  hashtablelist,      /**< pointer to hash list */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void*                 element             /**< element to append to the list */
   )
{
   SCIP_HASHTABLELIST* newlist;

   assert(hashtablelist != NULL);
   assert(blkmem != NULL);
   assert(element != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &newlist) );
   newlist->element = element;
   newlist->next = *hashtablelist;
   *hashtablelist = newlist;

   return SCIP_OKAY;
}

/** frees a hash list entry and all its successors */
static
void hashtablelistFree(
   SCIP_HASHTABLELIST**  hashtablelist,      /**< pointer to hash list to free */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_HASHTABLELIST* list;
   SCIP_HASHTABLELIST* nextlist;

   assert(hashtablelist != NULL);
   assert(blkmem != NULL);
   
   list = *hashtablelist;
   while( list != NULL )
   {
      nextlist = list->next;
      BMSfreeBlockMemory(blkmem, &list);
      list = nextlist;
   }

   *hashtablelist = NULL;
}

/** finds hash list entry pointing to element with given key in the hash list, returns NULL if not found */
static
SCIP_HASHTABLELIST* hashtablelistFind(
   SCIP_HASHTABLELIST*   hashtablelist,      /**< hash list */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr,            /**< user pointer */
   unsigned int          keyval,             /**< hash value of key */
   void*                 key                 /**< key to retrieve */
   )
{
   unsigned int currentkeyval;
   void* currentkey;

   assert(hashkeyeq != NULL);
   assert(key != NULL);

   while( hashtablelist != NULL )
   {
      currentkey = hashgetkey(userptr, hashtablelist->element);
      currentkeyval = hashkeyval(userptr, currentkey);
      if( currentkeyval == keyval && hashkeyeq(userptr, currentkey, key) )
         return hashtablelist;
      
      hashtablelist = hashtablelist->next;
   }

   return NULL;
}

/** retrieves element with given key from the hash list, or NULL */
static
void* hashtablelistRetrieve(
   SCIP_HASHTABLELIST*   hashtablelist,      /**< hash list */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr,            /**< user pointer */
   unsigned int          keyval,             /**< hash value of key */
   void*                 key                 /**< key to retrieve */
   )
{
   SCIP_HASHTABLELIST* h;

   /* find hash list entry */
   h = hashtablelistFind(hashtablelist, hashgetkey, hashkeyeq, hashkeyval, userptr, keyval, key);

   /* return element */
   if( h != NULL )
   {
#ifndef NDEBUG
      if( hashtablelistFind(h->next, hashgetkey, hashkeyeq, hashkeyval, userptr, keyval, key) != NULL )
      {
         SCIPwarningMessage("hashkey with same value exists multiple times (e.g. duplicate constraint/variable names), so the return value is maybe not correct\n");
      }
#endif
      
      return h->element;
   }
   else
      return NULL;
}


/** retrieves element with given key from the hash list, or NULL
 * returns pointer to hash table list entry */
static
void* hashtablelistRetrieveNext(
   SCIP_HASHTABLELIST**  hashtablelist,      /**< on input: hash list to search; on exit: hash list entry corresponding to element after retrieved one, or NULL */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr,            /**< user pointer */
   unsigned int          keyval,             /**< hash value of key */
   void*                 key                 /**< key to retrieve */
   )
{
   SCIP_HASHTABLELIST* h;

   assert(hashtablelist != NULL);
   
   /* find hash list entry */
   h = hashtablelistFind(*hashtablelist, hashgetkey, hashkeyeq, hashkeyval, userptr, keyval, key);

   /* return element */
   if( h != NULL )
   {
      *hashtablelist = h->next;
      
      return h->element;
   }
   
   *hashtablelist = NULL;
   
   return NULL;
}

/** removes element from the hash list */
static
SCIP_RETCODE hashtablelistRemove(
   SCIP_HASHTABLELIST**  hashtablelist,      /**< pointer to hash list */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void*                 element             /**< element to remove from the list */
   )
{
   SCIP_HASHTABLELIST* nextlist;

   assert(hashtablelist != NULL);
   assert(blkmem != NULL);
   assert(element != NULL);

   while( *hashtablelist != NULL && (*hashtablelist)->element != element )
      hashtablelist = &(*hashtablelist)->next;

   if( *hashtablelist != NULL )
   {
      nextlist = (*hashtablelist)->next;
      BMSfreeBlockMemory(blkmem, hashtablelist);
      *hashtablelist = nextlist;
   }

   return SCIP_OKAY;
}

/** creates a hash table */
SCIP_RETCODE SCIPhashtableCreate(
   SCIP_HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash table entries */
   int                   tablesize,          /**< size of the hash table */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr             /**< user pointer */
   )
{
   assert(hashtable != NULL);
   assert(tablesize > 0);
   assert(hashgetkey != NULL);
   assert(hashkeyeq != NULL);
   assert(hashkeyval != NULL);

   SCIP_ALLOC( BMSallocMemory(hashtable) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*hashtable)->lists, tablesize) );
   (*hashtable)->blkmem = blkmem;
   (*hashtable)->nlists = tablesize;
   (*hashtable)->hashgetkey = hashgetkey;
   (*hashtable)->hashkeyeq = hashkeyeq;
   (*hashtable)->hashkeyval = hashkeyval;
   (*hashtable)->userptr = userptr;

   BMSclearMemoryArray((*hashtable)->lists, tablesize);

   return SCIP_OKAY;
}

/** frees the hash table */
void SCIPhashtableFree(
   SCIP_HASHTABLE**      hashtable           /**< pointer to the hash table */
   )
{
   int i;
   SCIP_HASHTABLE* table;
   BMS_BLKMEM* blkmem;
   SCIP_HASHTABLELIST** lists;

   assert(hashtable != NULL);
   assert(*hashtable != NULL);

   table = (*hashtable);
   blkmem = table->blkmem;
   lists = table->lists;

   /* free hash lists */
   for( i = table->nlists - 1; i >= 0; --i )
      hashtablelistFree(&lists[i], blkmem);

   /* free main hast table data structure */
   BMSfreeMemoryArray(&table->lists);
   BMSfreeMemory(hashtable);
}

/** inserts element in hash table (multiple inserts of same element possible) */
SCIP_RETCODE SCIPhashtableInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   )
{
   void* key;
   unsigned int keyval;
   unsigned int hashval;

   assert(hashtable != NULL);
   assert(hashtable->lists != NULL);
   assert(hashtable->nlists > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(element != NULL);

   /* get the hash key and its hash value */
   key = hashtable->hashgetkey(hashtable->userptr, element);
   keyval = hashtable->hashkeyval(hashtable->userptr, key);
   hashval = keyval % hashtable->nlists; /*lint !e573*/

   /* append element to the list at the hash position */
   SCIP_CALL( hashtablelistAppend(&hashtable->lists[hashval], hashtable->blkmem, element) );
   
   return SCIP_OKAY;
}

/** inserts element in hash table (multiple insertion of same element is checked and results in an error) */
SCIP_RETCODE SCIPhashtableSafeInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   )
{
   assert(hashtable != NULL);
   assert(hashtable->hashgetkey != NULL);

   /* check, if key is already existing */
   if( SCIPhashtableRetrieve(hashtable, hashtable->hashgetkey(hashtable->userptr, element)) != NULL )
      return SCIP_KEYALREADYEXISTING;

   /* insert element in hash table */
   SCIP_CALL( SCIPhashtableInsert(hashtable, element) );
   
   return SCIP_OKAY;
}

/** retrieve element with key from hash table, returns NULL if not existing */
void* SCIPhashtableRetrieve(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 key                 /**< key to retrieve */
   )
{
   unsigned int keyval;
   unsigned int hashval;

   assert(hashtable != NULL);
   assert(hashtable->lists != NULL);
   assert(hashtable->nlists > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(key != NULL);

   /* get the hash value of the key */
   keyval = hashtable->hashkeyval(hashtable->userptr, key);
   hashval = keyval % hashtable->nlists; /*lint !e573*/

   return hashtablelistRetrieve(hashtable->lists[hashval], hashtable->hashgetkey, hashtable->hashkeyeq, 
      hashtable->hashkeyval, hashtable->userptr, keyval, key);
}

/** retrieve element with key from hash table, returns NULL if not existing
 * can be used to retrieve all entries with the same key (one-by-one) */
void* SCIPhashtableRetrieveNext(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   SCIP_HASHTABLELIST**  hashtablelist,      /**< input: entry in hash table list from which to start searching, or NULL; output: entry in hash table list corresponding to element after retrieved one, or NULL */
   void*                 key                 /**< key to retrieve */
   )
{
   unsigned int keyval;

   assert(hashtable != NULL);
   assert(hashtable->lists != NULL);
   assert(hashtable->nlists > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(hashtablelist != NULL);
   assert(key != NULL);

   keyval = hashtable->hashkeyval(hashtable->userptr, key);

   if( *hashtablelist == NULL )
   {
      unsigned int hashval;
      
      /* get the hash value of the key */
      hashval = keyval % hashtable->nlists; /*lint !e573*/
      
      *hashtablelist = hashtable->lists[hashval];
   }

   return hashtablelistRetrieveNext(hashtablelist, hashtable->hashgetkey, hashtable->hashkeyeq, 
      hashtable->hashkeyval, hashtable->userptr, keyval, key);
}

/** returns whether the given element exists in the table */
SCIP_Bool SCIPhashtableExists(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to search in the table */
   )
{
   void* key;
   unsigned int keyval;
   unsigned int hashval;

   assert(hashtable != NULL);
   assert(hashtable->lists != NULL);
   assert(hashtable->nlists > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(element != NULL);

   /* get the hash key and its hash value */
   key = hashtable->hashgetkey(hashtable->userptr, element);
   keyval = hashtable->hashkeyval(hashtable->userptr, key);
   hashval = keyval % hashtable->nlists; /*lint !e573*/

   return (hashtablelistFind(hashtable->lists[hashval], hashtable->hashgetkey, hashtable->hashkeyeq,
         hashtable->hashkeyval, hashtable->userptr, keyval, key) != NULL);
}

/** removes element from the hash table, if it exists */
SCIP_RETCODE SCIPhashtableRemove(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to remove from the table */
   )
{
   void* key;
   unsigned int keyval;
   unsigned int hashval;

   assert(hashtable != NULL);
   assert(hashtable->lists != NULL);
   assert(hashtable->nlists > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(element != NULL);

   /* get the hash key and its hash value */
   key = hashtable->hashgetkey(hashtable->userptr, element);
   keyval = hashtable->hashkeyval(hashtable->userptr, key);
   hashval = keyval % hashtable->nlists; /*lint !e573*/

   /* remove element from the list at the hash position */
   SCIP_CALL( hashtablelistRemove(&hashtable->lists[hashval], hashtable->blkmem, element) );
   
   return SCIP_OKAY;
}

/** prints statistics about hash table usage */
void SCIPhashtablePrintStatistics(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   )
{
   SCIP_HASHTABLELIST* hashtablelist;
   int usedslots;
   int maxslotsize;
   int sumslotsize;
   int slotsize;
   int i;

   assert(hashtable != NULL);

   usedslots = 0;
   maxslotsize = 0;
   sumslotsize = 0;
   for( i = 0; i < hashtable->nlists; ++i )
   {
      hashtablelist = hashtable->lists[i];
      if( hashtablelist != NULL )
      {
         usedslots++;
         slotsize = 0;
         while( hashtablelist != NULL )
         {
            slotsize++;
            hashtablelist = hashtablelist->next;
         }
         maxslotsize = MAX(maxslotsize, slotsize);
         sumslotsize += slotsize;
      }
   }

   SCIPmessagePrintInfo("%d hash entries, used %d/%d slots (%.1f%%)",
      sumslotsize, usedslots, hashtable->nlists, 100.0*(SCIP_Real)usedslots/(SCIP_Real)(hashtable->nlists));
   if( usedslots > 0 )
      SCIPmessagePrintInfo(", avg. %.1f entries/used slot, max. %d entries in slot",
         (SCIP_Real)sumslotsize/(SCIP_Real)usedslots, maxslotsize);
   SCIPmessagePrintInfo("\n");
}


/** returns TRUE iff both keys (i.e. strings) are equal */
SCIP_DECL_HASHKEYEQ(SCIPhashKeyEqString)
{  /*lint --e{715}*/
   const char* string1 = (const char*)key1;
   const char* string2 = (const char*)key2;

   return (strcmp(string1, string2) == 0);
}

/** returns the hash value of the key (i.e. string) */
SCIP_DECL_HASHKEYVAL(SCIPhashKeyValString)
{  /*lint --e{715}*/
   const char* str;
   unsigned int sum;

   str = (const char*)key;
   sum = 0;
   while( *str != '\0' )
   {
      sum *= 31;
      sum += (unsigned int)(*str); /*lint !e571*/
      str++;
   }

   return sum;
}




/*
 * Hash Map
 */

/** appends origin->image pair to the hash list */
static
SCIP_RETCODE hashmaplistAppend(
   SCIP_HASHMAPLIST**    hashmaplist,        /**< pointer to hash list */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL */
   void*                 origin,             /**< origin of the mapping origin -> image */
   void*                 image               /**< image of the mapping origin -> image */
   )
{
   SCIP_HASHMAPLIST* newlist;

   assert(hashmaplist != NULL);
   assert(origin != NULL);

   if( blkmem != NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &newlist) );
   }
   else
   {
      SCIP_ALLOC( BMSallocMemory(&newlist) );
   }

   newlist->origin = origin;
   newlist->image = image;
   newlist->next = *hashmaplist;
   *hashmaplist = newlist;

   return SCIP_OKAY;
}

/** frees a hash list entry and all its successors */
static
void hashmaplistFree(
   SCIP_HASHMAPLIST**    hashmaplist,        /**< pointer to hash list to free */
   BMS_BLKMEM*           blkmem              /**< block memory, or NULL */
   )
{
   SCIP_HASHMAPLIST* list;
   SCIP_HASHMAPLIST* nextlist;

   assert(hashmaplist != NULL);
   
   list = *hashmaplist;
   while( list != NULL )
   {
      nextlist = list->next;

      if( blkmem != NULL )
      {
         BMSfreeBlockMemory(blkmem, &list);
      }
      else
      {
         BMSfreeMemory(&list);
      }

      list = nextlist;
   }

   *hashmaplist = NULL;
}

/** finds hash list entry pointing to given origin in the hash list, returns NULL if not found */
static
SCIP_HASHMAPLIST* hashmaplistFind(
   SCIP_HASHMAPLIST*     hashmaplist,        /**< hash list */
   void*                 origin              /**< origin to find */
   )
{
   assert(origin != NULL);

   while( hashmaplist != NULL )
   {
      if( hashmaplist->origin == origin )
         return hashmaplist;
      hashmaplist = hashmaplist->next;
   }

   return NULL;
}

/** retrieves image of given origin from the hash list, or NULL */
static
void* hashmaplistGetImage(
   SCIP_HASHMAPLIST*     hashmaplist,        /**< hash list */
   void*                 origin              /**< origin to retrieve image for */
   )
{
   SCIP_HASHMAPLIST* h;

   /* find hash list entry */
   h = hashmaplistFind(hashmaplist, origin);

   /* return image */
   if( h != NULL )
      return h->image;
   else
      return NULL;
}

/** sets image for given origin in the hash list, either by modifying existing origin->image pair or by appending a
 *  new origin->image pair
 */
static
SCIP_RETCODE hashmaplistSetImage(
   SCIP_HASHMAPLIST**    hashmaplist,        /**< pointer to hash list */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   )
{
   SCIP_HASHMAPLIST* h;

   /* find hash list entry */
   h = hashmaplistFind(*hashmaplist, origin);

   /* set image or add origin->image pair */
   if( h != NULL )
      h->image = image;
   else
   {
      SCIP_CALL( hashmaplistAppend(hashmaplist, blkmem, origin, image) );
   }

   return SCIP_OKAY;
}

/** removes origin->image pair from the hash list */
static
SCIP_RETCODE hashmaplistRemove(
   SCIP_HASHMAPLIST**    hashmaplist,        /**< pointer to hash list */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL */
   void*                 origin              /**< origin to remove from the list */
   )
{
   SCIP_HASHMAPLIST* nextlist;

   assert(hashmaplist != NULL);
   assert(origin != NULL);

   while( *hashmaplist != NULL && (*hashmaplist)->origin != origin )
   {
      hashmaplist = &(*hashmaplist)->next;
   }
   if( *hashmaplist != NULL )
   {
      nextlist = (*hashmaplist)->next;

      if( blkmem != NULL )
      {
         BMSfreeBlockMemory(blkmem, hashmaplist);
      }
      else
      {
         BMSfreeMemory(hashmaplist);
      }

      *hashmaplist = nextlist;
   }

   return SCIP_OKAY;
}


/** creates a hash map mapping pointers to pointers 
 *
 * @note if possible always use a blkmem pointer instead of NULL, otherwise it could slow down the map
 */
SCIP_RETCODE SCIPhashmapCreate(
   SCIP_HASHMAP**        hashmap,            /**< pointer to store the created hash map */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash map entries, or NULL */
   int                   mapsize             /**< size of the hash map */
   )
{
   int i;

   assert(hashmap != NULL);
   assert(mapsize > 0);

   SCIP_ALLOC( BMSallocMemory(hashmap) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*hashmap)->lists, mapsize) );
   (*hashmap)->blkmem = blkmem;
   (*hashmap)->nlists = mapsize;

   /* initialize hash lists */
   for( i = 0; i < mapsize; ++i )
      (*hashmap)->lists[i] = NULL;

   return SCIP_OKAY;
}

/** frees the hash map */
void SCIPhashmapFree(
   SCIP_HASHMAP**        hashmap             /**< pointer to the hash map */
   )
{
   int i;

   assert(hashmap != NULL);
   assert(*hashmap != NULL);

   /* free hash lists */
   for( i = 0; i < (*hashmap)->nlists; ++i )
      hashmaplistFree(&(*hashmap)->lists[i], (*hashmap)->blkmem);

   /* free main hash map data structure */
   BMSfreeMemoryArray(&(*hashmap)->lists);
   BMSfreeMemory(hashmap);
}

/** inserts new origin->image pair in hash map (must not be called for already existing origins!) */
SCIP_RETCODE SCIPhashmapInsert(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)((size_t)origin % (unsigned int)hashmap->nlists);

   /* append origin->image pair to the list at the hash position */
   SCIP_CALL( hashmaplistAppend(&hashmap->lists[hashval], hashmap->blkmem, origin, image) );
   
   return SCIP_OKAY;
}

/** retrieves image of given origin from the hash map, or NULL if no image exists */
void* SCIPhashmapGetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to retrieve image for */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)((size_t)origin % (unsigned int)hashmap->nlists);

   /* get image for origin from hash list */
   return hashmaplistGetImage(hashmap->lists[hashval], origin);
}

/** sets image for given origin in the hash map, either by modifying existing origin->image pair or by appending a
 *  new origin->image pair
 */
SCIP_RETCODE SCIPhashmapSetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)((size_t)origin % (unsigned int)hashmap->nlists);

   /* set image for origin in hash list */
   SCIP_CALL( hashmaplistSetImage(&hashmap->lists[hashval], hashmap->blkmem, origin, image) );
   
   return SCIP_OKAY;
}

/** checks whether an image to the given origin exists in the hash map */
SCIP_Bool SCIPhashmapExists(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to search for */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)((size_t)origin % (unsigned int)hashmap->nlists);

   return (hashmaplistFind(hashmap->lists[hashval], origin) != NULL);
}

/** removes origin->image pair from the hash map, if it exists */
SCIP_RETCODE SCIPhashmapRemove(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to remove from the list */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)((size_t)origin % (unsigned int)hashmap->nlists);

   /* remove element from the list at the hash position */
   SCIP_CALL( hashmaplistRemove(&hashmap->lists[hashval], hashmap->blkmem, origin) );
   
   return SCIP_OKAY;
}

/** prints statistics about hash map usage */
void SCIPhashmapPrintStatistics(
   SCIP_HASHMAP*         hashmap             /**< hash map */
   )
{
   SCIP_HASHMAPLIST* hashmaplist;
   int usedslots;
   int maxslotsize;
   int sumslotsize;
   int slotsize;
   int i;

   assert(hashmap != NULL);

   usedslots = 0;
   maxslotsize = 0;
   sumslotsize = 0;
   for( i = 0; i < hashmap->nlists; ++i )
   {
      hashmaplist = hashmap->lists[i];
      if( hashmaplist != NULL )
      {
         usedslots++;
         slotsize = 0;
         while( hashmaplist != NULL )
         {
            slotsize++;
            hashmaplist = hashmaplist->next;
         }
         maxslotsize = MAX(maxslotsize, slotsize);
         sumslotsize += slotsize;
      }
   }

   SCIPmessagePrintInfo("%d hash entries, used %d/%d slots (%.1f%%)",
      sumslotsize, usedslots, hashmap->nlists, 100.0*(SCIP_Real)usedslots/(SCIP_Real)(hashmap->nlists));
   if( usedslots > 0 )
      SCIPmessagePrintInfo(", avg. %.1f entries/used slot, max. %d entries in slot", 
         (SCIP_Real)sumslotsize/(SCIP_Real)usedslots, maxslotsize);
   SCIPmessagePrintInfo("\n");
}

/** indicates whether a hash map has no entries */
SCIP_Bool SCIPhashmapIsEmpty(
   SCIP_HASHMAP*      hashmap          /**< hash map */
)
{
   int i;
   assert(hashmap != NULL);
   
   for( i = 0; i < hashmap->nlists; ++i )
      if( hashmap->lists[i] )
         return FALSE;
   
   return TRUE;
}

/** gives the number of entries in a hash map */ 
int SCIPhashmapGetNEntries(
   SCIP_HASHMAP*      hashmap          /**< hash map */
)
{
   int count = 0;
   int i;
   assert(hashmap != NULL);
   
   for( i = 0; i < hashmap->nlists; ++i )
      count += SCIPhashmapListGetNEntries(hashmap->lists[i]);

   return count;
}

/** gives the number of lists (buckets) in a hash map */ 
int SCIPhashmapGetNLists(
   SCIP_HASHMAP*      hashmap          /**< hash map */
)
{
   assert(hashmap != NULL);
   
   return hashmap->nlists;
}

/** gives a specific list (bucket) in a hash map */
SCIP_HASHMAPLIST* SCIPhashmapGetList(
   SCIP_HASHMAP*     hashmap,          /**< hash map */
   int               listindex         /**< index of hash map list */
)
{
   assert(hashmap != NULL);
   assert(listindex >= 0);
   assert(listindex < hashmap->nlists);
   
   return hashmap->lists[listindex];
}

/** gives the number of entries in a list of a hash map */ 
int SCIPhashmapListGetNEntries(
   SCIP_HASHMAPLIST* hashmaplist       /**< hash map list, can be NULL */
)
{
   int count = 0;
   
   for( ; hashmaplist; hashmaplist = hashmaplist->next )
      ++count;
   
   return count;
}

/** retrieves origin of given entry in a hash map */ 
void* SCIPhashmapListGetOrigin(
   SCIP_HASHMAPLIST* hashmaplist       /**< hash map list */
)
{
   assert(hashmaplist != NULL);
   
   return hashmaplist->origin;
}

/** retrieves image of given entry in a hash map */ 
void* SCIPhashmapListGetImage(
   SCIP_HASHMAPLIST* hashmaplist       /**< hash map list */
)
{
   assert(hashmaplist != NULL);
   
   return hashmaplist->image;
}

/** retrieves next entry from given entry in a hash map list, or NULL if at end of list. */ 
SCIP_HASHMAPLIST* SCIPhashmapListGetNext(
   SCIP_HASHMAPLIST* hashmaplist       /**< hash map list */
)
{
   assert(hashmaplist != NULL);
   
   return hashmaplist->next;
}

/** removes all entries in a hash map. */ 
SCIP_RETCODE SCIPhashmapRemoveAll(
   SCIP_HASHMAP*     hashmap           /**< hash map */
)
{
   int listidx;

   assert(hashmap != NULL);
   
   /* free hash lists */
   for( listidx = 0; listidx < hashmap->nlists; ++listidx )
      hashmaplistFree(&hashmap->lists[listidx], hashmap->blkmem);

   return SCIP_OKAY;
}



/*
 * Dynamic Arrays
 */

/** creates a dynamic array of real values */
SCIP_RETCODE SCIPrealarrayCreate(
   SCIP_REALARRAY**      realarray,          /**< pointer to store the real array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(realarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, realarray) );
   (*realarray)->blkmem = blkmem;
   (*realarray)->vals = NULL;
   (*realarray)->valssize = 0;
   (*realarray)->firstidx = -1;
   (*realarray)->minusedidx = INT_MAX;
   (*realarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of real values */
SCIP_RETCODE SCIPrealarrayCopy(
   SCIP_REALARRAY**      realarray,          /**< pointer to store the copied real array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REALARRAY*       sourcerealarray     /**< dynamic real array to copy */
   )
{
   assert(realarray != NULL);
   assert(sourcerealarray != NULL);

   SCIP_CALL( SCIPrealarrayCreate(realarray, blkmem) );
   if( sourcerealarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*realarray)->vals, sourcerealarray->vals,
                     sourcerealarray->valssize) );
   }
   (*realarray)->valssize = sourcerealarray->valssize;
   (*realarray)->firstidx = sourcerealarray->firstidx;
   (*realarray)->minusedidx = sourcerealarray->minusedidx;
   (*realarray)->maxusedidx = sourcerealarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of real values */
SCIP_RETCODE SCIPrealarrayFree(
   SCIP_REALARRAY**      realarray           /**< pointer to the real array */
   )
{
   assert(realarray != NULL);
   assert(*realarray != NULL);

   BMSfreeBlockMemoryArrayNull((*realarray)->blkmem, &(*realarray)->vals, (*realarray)->valssize);
   BMSfreeBlockMemory((*realarray)->blkmem, realarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPrealarrayExtend(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(realarray != NULL);
   assert(realarray->minusedidx == INT_MAX || realarray->firstidx >= 0);
   assert(realarray->maxusedidx == INT_MIN || realarray->firstidx >= 0);
   assert(realarray->minusedidx == INT_MAX || realarray->minusedidx >= realarray->firstidx);
   assert(realarray->maxusedidx == INT_MIN || realarray->maxusedidx < realarray->firstidx + realarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);
   
   minidx = MIN(minidx, realarray->minusedidx);
   maxidx = MAX(maxidx, realarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   SCIPdebugMessage("extending realarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n", 
      (void*)realarray, realarray->firstidx, realarray->valssize, realarray->minusedidx, realarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > realarray->valssize )
   {
      SCIP_Real* newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = SCIPsetCalcMemGrowSize(set, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(realarray->blkmem, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( realarray->firstidx != -1 )
      {
         for( i = 0; i < realarray->minusedidx - newfirstidx; ++i )
            newvals[i] = 0.0;

         /* check for possible overflow or negative value */
         assert(realarray->maxusedidx - realarray->minusedidx + 1 > 0);

         BMScopyMemoryArray(&newvals[realarray->minusedidx - newfirstidx],
            &(realarray->vals[realarray->minusedidx - realarray->firstidx]),
            realarray->maxusedidx - realarray->minusedidx + 1); /*lint !e866*/
         for( i = realarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = 0.0;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = 0.0;
      }

      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(realarray->blkmem, &realarray->vals, realarray->valssize);
      realarray->vals = newvals;
      realarray->valssize = newvalssize;
      realarray->firstidx = newfirstidx;
   }
   else if( realarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = realarray->valssize - nused;
      assert(nfree >= 0);
      realarray->firstidx = minidx - nfree/2;
      assert(realarray->firstidx <= minidx);
      assert(maxidx < realarray->firstidx + realarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < realarray->valssize; ++i )
         assert(realarray->vals[i] == 0.0);
#endif
   }
   else if( minidx < realarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = realarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + realarray->valssize);
      
      if( realarray->minusedidx <= realarray->maxusedidx )
      {
         int shift;

         assert(realarray->firstidx <= realarray->minusedidx);
         assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);

         /* shift used part of array to the right */
         shift = realarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = realarray->maxusedidx - realarray->firstidx; i >= realarray->minusedidx - realarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < realarray->valssize);
            realarray->vals[i + shift] = realarray->vals[i];
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            realarray->vals[realarray->minusedidx - realarray->firstidx + i] = 0.0;
      }
      realarray->firstidx = newfirstidx;
   }
   else if( maxidx >= realarray->firstidx + realarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = realarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + realarray->valssize);
      
      if( realarray->minusedidx <= realarray->maxusedidx )
      {
         int shift;

         assert(realarray->firstidx <= realarray->minusedidx);
         assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - realarray->firstidx;
         assert(shift > 0);
         for( i = realarray->minusedidx - realarray->firstidx; i <= realarray->maxusedidx - realarray->firstidx; ++i )
         {
            assert(0 <= i - shift && i - shift < realarray->valssize);
            realarray->vals[i - shift] = realarray->vals[i];
         }
         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            realarray->vals[realarray->maxusedidx - realarray->firstidx - i] = 0.0;
      }
      realarray->firstidx = newfirstidx;
   }

   assert(minidx >= realarray->firstidx);
   assert(maxidx < realarray->firstidx + realarray->valssize);

   return SCIP_OKAY;
}

/** clears a dynamic real array */
SCIP_RETCODE SCIPrealarrayClear(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(realarray != NULL);

   SCIPdebugMessage("clearing realarray %p (firstidx=%d, size=%d, range=[%d,%d])\n", 
      (void*)realarray, realarray->firstidx, realarray->valssize, realarray->minusedidx, realarray->maxusedidx);

   if( realarray->minusedidx <= realarray->maxusedidx )
   {
      int i;
   
      assert(realarray->firstidx <= realarray->minusedidx);
      assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);
      assert(realarray->firstidx != -1);
      assert(realarray->valssize > 0);

      /* clear the used part of array */
      for( i = realarray->minusedidx - realarray->firstidx; i <= realarray->maxusedidx - realarray->firstidx; ++i )
         realarray->vals[i] = 0.0;

      /* mark the array cleared */
      realarray->minusedidx = INT_MAX;
      realarray->maxusedidx = INT_MIN;
   }
   assert(realarray->minusedidx == INT_MAX);
   assert(realarray->maxusedidx == INT_MIN);

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
SCIP_Real SCIPrealarrayGetVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(realarray != NULL);
   assert(idx >= 0);
   
   if( idx < realarray->minusedidx || idx > realarray->maxusedidx )
      return 0.0;
   else
   {
      assert(realarray->vals != NULL);
      assert(idx - realarray->firstidx >= 0);
      assert(idx - realarray->firstidx < realarray->valssize);

      return realarray->vals[idx - realarray->firstidx];
   }
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPrealarraySetVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   idx,                /**< array index to set value for */
   SCIP_Real             val                 /**< value to set array index to */
   )
{
   assert(realarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting realarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %g\n", 
      (void*)realarray, realarray->firstidx, realarray->valssize, realarray->minusedidx, realarray->maxusedidx, idx, val);

   if( !SCIPsetIsZero(set, val) )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPrealarrayExtend(realarray, set, idx, idx) );
      assert(idx >= realarray->firstidx);
      assert(idx < realarray->firstidx + realarray->valssize);
      
      /* set the array value of the index */
      realarray->vals[idx - realarray->firstidx] = val;

      /* update min/maxusedidx */
      realarray->minusedidx = MIN(realarray->minusedidx, idx);
      realarray->maxusedidx = MAX(realarray->maxusedidx, idx);
   }
   else if( idx >= realarray->firstidx && idx < realarray->firstidx + realarray->valssize )
   {
      /* set the array value of the index to zero */
      realarray->vals[idx - realarray->firstidx] = 0.0;
      
      /* check, if we can tighten the min/maxusedidx */
      if( idx == realarray->minusedidx )
      {
         assert(realarray->maxusedidx >= 0);
         assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);
         do
         {
            realarray->minusedidx++;
         }
         while( realarray->minusedidx <= realarray->maxusedidx
            && SCIPsetIsZero(set, realarray->vals[realarray->minusedidx - realarray->firstidx]) );
         if( realarray->minusedidx > realarray->maxusedidx )
         {
            realarray->minusedidx = INT_MAX;
            realarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == realarray->maxusedidx )
      {
         assert(realarray->minusedidx >= 0);
         assert(realarray->minusedidx < realarray->maxusedidx);
         assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);
         do
         {
            realarray->maxusedidx--;
            assert(realarray->minusedidx <= realarray->maxusedidx);
         }
         while( SCIPsetIsZero(set, realarray->vals[realarray->maxusedidx - realarray->firstidx]) );
      }      
   }

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
SCIP_RETCODE SCIPrealarrayIncVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   idx,                /**< array index to increase value for */
   SCIP_Real             incval              /**< value to increase array index */
   )
{
   SCIP_Real oldval;

   oldval = SCIPrealarrayGetVal(realarray, idx);
   if( oldval != SCIP_INVALID ) /*lint !e777*/
      return SCIPrealarraySetVal(realarray, set, idx, oldval + incval);
   else
      return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPrealarrayGetMinIdx(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(realarray != NULL);

   return realarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPrealarrayGetMaxIdx(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(realarray != NULL);

   return realarray->maxusedidx;
}

/** creates a dynamic array of int values */
SCIP_RETCODE SCIPintarrayCreate(
   SCIP_INTARRAY**       intarray,           /**< pointer to store the int array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(intarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, intarray) );
   (*intarray)->blkmem = blkmem;
   (*intarray)->vals = NULL;
   (*intarray)->valssize = 0;
   (*intarray)->firstidx = -1;
   (*intarray)->minusedidx = INT_MAX;
   (*intarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of int values */
SCIP_RETCODE SCIPintarrayCopy(
   SCIP_INTARRAY**       intarray,           /**< pointer to store the copied int array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_INTARRAY*        sourceintarray      /**< dynamic int array to copy */
   )
{
   assert(intarray != NULL);
   assert(sourceintarray != NULL);

   SCIP_CALL( SCIPintarrayCreate(intarray, blkmem) );
   if( sourceintarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*intarray)->vals, sourceintarray->vals, sourceintarray->valssize) );
   }
   (*intarray)->valssize = sourceintarray->valssize;
   (*intarray)->firstidx = sourceintarray->firstidx;
   (*intarray)->minusedidx = sourceintarray->minusedidx;
   (*intarray)->maxusedidx = sourceintarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of int values */
SCIP_RETCODE SCIPintarrayFree(
   SCIP_INTARRAY**       intarray            /**< pointer to the int array */
   )
{
   assert(intarray != NULL);
   assert(*intarray != NULL);

   BMSfreeBlockMemoryArrayNull((*intarray)->blkmem, &(*intarray)->vals, (*intarray)->valssize);
   BMSfreeBlockMemory((*intarray)->blkmem, intarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPintarrayExtend(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(intarray != NULL);
   assert(intarray->minusedidx == INT_MAX || intarray->firstidx >= 0);
   assert(intarray->maxusedidx == INT_MIN || intarray->firstidx >= 0);
   assert(intarray->minusedidx == INT_MAX || intarray->minusedidx >= intarray->firstidx);
   assert(intarray->maxusedidx == INT_MIN || intarray->maxusedidx < intarray->firstidx + intarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);
   
   minidx = MIN(minidx, intarray->minusedidx);
   maxidx = MAX(maxidx, intarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   SCIPdebugMessage("extending intarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n", 
      (void*)intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > intarray->valssize )
   {
      int* newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = SCIPsetCalcMemGrowSize(set, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(intarray->blkmem, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( intarray->firstidx != -1 )
      {
         for( i = 0; i < intarray->minusedidx - newfirstidx; ++i )
            newvals[i] = 0;

         /* check for possible overflow or negative value */
         assert(intarray->maxusedidx - intarray->minusedidx + 1 > 0);

         BMScopyMemoryArray(&newvals[intarray->minusedidx - newfirstidx],
            &intarray->vals[intarray->minusedidx - intarray->firstidx],
            intarray->maxusedidx - intarray->minusedidx + 1); /*lint !e866*/
         for( i = intarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = 0;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = 0;
      }

      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(intarray->blkmem, &intarray->vals, intarray->valssize);
      intarray->vals = newvals;
      intarray->valssize = newvalssize;
      intarray->firstidx = newfirstidx;
   }
   else if( intarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = intarray->valssize - nused;
      assert(nfree >= 0);
      intarray->firstidx = minidx - nfree/2;
      assert(intarray->firstidx <= minidx);
      assert(maxidx < intarray->firstidx + intarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < intarray->valssize; ++i )
         assert(intarray->vals[i] == 0);
#endif
   }
   else if( minidx < intarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = intarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + intarray->valssize);
      
      if( intarray->minusedidx <= intarray->maxusedidx )
      {
         int shift;

         assert(intarray->firstidx <= intarray->minusedidx);
         assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);

         /* shift used part of array to the right */
         shift = intarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = intarray->maxusedidx - intarray->firstidx; i >= intarray->minusedidx - intarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < intarray->valssize);
            intarray->vals[i + shift] = intarray->vals[i];
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            intarray->vals[intarray->minusedidx - intarray->firstidx + i] = 0;
      }
      intarray->firstidx = newfirstidx;
   }
   else if( maxidx >= intarray->firstidx + intarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = intarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + intarray->valssize);
      
      if( intarray->minusedidx <= intarray->maxusedidx )
      {
         int shift;

         assert(intarray->firstidx <= intarray->minusedidx);
         assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - intarray->firstidx;
         assert(shift > 0);
         for( i = intarray->minusedidx - intarray->firstidx; i <= intarray->maxusedidx - intarray->firstidx; ++i )
         {
            assert(0 <= i - shift && i - shift < intarray->valssize);
            intarray->vals[i - shift] = intarray->vals[i];
         }
         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            intarray->vals[intarray->maxusedidx - intarray->firstidx - i] = 0;
      }
      intarray->firstidx = newfirstidx;
   }

   assert(minidx >= intarray->firstidx);
   assert(maxidx < intarray->firstidx + intarray->valssize);

   return SCIP_OKAY;
}

/** clears a dynamic int array */
SCIP_RETCODE SCIPintarrayClear(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(intarray != NULL);

   SCIPdebugMessage("clearing intarray %p (firstidx=%d, size=%d, range=[%d,%d])\n", 
      (void*)intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx);

   if( intarray->minusedidx <= intarray->maxusedidx )
   {
      int i;
   
      assert(intarray->firstidx <= intarray->minusedidx);
      assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);
      assert(intarray->firstidx != -1);
      assert(intarray->valssize > 0);

      /* clear the used part of array */
      for( i = intarray->minusedidx - intarray->firstidx; i <= intarray->maxusedidx - intarray->firstidx; ++i )
         intarray->vals[i] = 0;

      /* mark the array cleared */
      intarray->minusedidx = INT_MAX;
      intarray->maxusedidx = INT_MIN;
   }
   assert(intarray->minusedidx == INT_MAX);
   assert(intarray->maxusedidx == INT_MIN);

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
int SCIPintarrayGetVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(intarray != NULL);
   assert(idx >= 0);
   
   if( idx < intarray->minusedidx || idx > intarray->maxusedidx )
      return 0;
   else
   {
      assert(intarray->vals != NULL);
      assert(idx - intarray->firstidx >= 0);
      assert(idx - intarray->firstidx < intarray->valssize);

      return intarray->vals[idx - intarray->firstidx];
   }
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPintarraySetVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   idx,                /**< array index to set value for */
   int                   val                 /**< value to set array index to */
   )
{
   assert(intarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting intarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %d\n", 
      (void*)intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx, idx, val);

   if( val != 0 )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPintarrayExtend(intarray, set, idx, idx) );
      assert(idx >= intarray->firstidx);
      assert(idx < intarray->firstidx + intarray->valssize);
      
      /* set the array value of the index */
      intarray->vals[idx - intarray->firstidx] = val;

      /* update min/maxusedidx */
      intarray->minusedidx = MIN(intarray->minusedidx, idx);
      intarray->maxusedidx = MAX(intarray->maxusedidx, idx);
   }
   else if( idx >= intarray->firstidx && idx < intarray->firstidx + intarray->valssize )
   {
      /* set the array value of the index to zero */
      intarray->vals[idx - intarray->firstidx] = 0;
      
      /* check, if we can tighten the min/maxusedidx */
      if( idx == intarray->minusedidx )
      {
         assert(intarray->maxusedidx >= 0);
         assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);
         do
         {
            intarray->minusedidx++;
         }
         while( intarray->minusedidx <= intarray->maxusedidx
            && intarray->vals[intarray->minusedidx - intarray->firstidx] == 0 );
         if( intarray->minusedidx > intarray->maxusedidx )
         {
            intarray->minusedidx = INT_MAX;
            intarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == intarray->maxusedidx )
      {
         assert(intarray->minusedidx >= 0);
         assert(intarray->minusedidx < intarray->maxusedidx);
         assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);
         do
         {
            intarray->maxusedidx--;
            assert(intarray->minusedidx <= intarray->maxusedidx);
         }
         while( intarray->vals[intarray->maxusedidx - intarray->firstidx] == 0 );
      }      
   }

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
SCIP_RETCODE SCIPintarrayIncVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   idx,                /**< array index to increase value for */
   int                   incval              /**< value to increase array index */
   )
{
   return SCIPintarraySetVal(intarray, set, idx, SCIPintarrayGetVal(intarray, idx) + incval);
}

/** returns the minimal index of all stored non-zero elements */
int SCIPintarrayGetMinIdx(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(intarray != NULL);

   return intarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPintarrayGetMaxIdx(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(intarray != NULL);

   return intarray->maxusedidx;
}


/** creates a dynamic array of bool values */
SCIP_RETCODE SCIPboolarrayCreate(
   SCIP_BOOLARRAY**      boolarray,          /**< pointer to store the bool array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(boolarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, boolarray) );
   (*boolarray)->blkmem = blkmem;
   (*boolarray)->vals = NULL;
   (*boolarray)->valssize = 0;
   (*boolarray)->firstidx = -1;
   (*boolarray)->minusedidx = INT_MAX;
   (*boolarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of bool values */
SCIP_RETCODE SCIPboolarrayCopy(
   SCIP_BOOLARRAY**      boolarray,          /**< pointer to store the copied bool array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_BOOLARRAY*       sourceboolarray     /**< dynamic bool array to copy */
   )
{
   assert(boolarray != NULL);
   assert(sourceboolarray != NULL);

   SCIP_CALL( SCIPboolarrayCreate(boolarray, blkmem) );
   if( sourceboolarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*boolarray)->vals, sourceboolarray->vals, 
                     sourceboolarray->valssize) );
   }
   (*boolarray)->valssize = sourceboolarray->valssize;
   (*boolarray)->firstidx = sourceboolarray->firstidx;
   (*boolarray)->minusedidx = sourceboolarray->minusedidx;
   (*boolarray)->maxusedidx = sourceboolarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of bool values */
SCIP_RETCODE SCIPboolarrayFree(
   SCIP_BOOLARRAY**      boolarray           /**< pointer to the bool array */
   )
{
   assert(boolarray != NULL);
   assert(*boolarray != NULL);

   BMSfreeBlockMemoryArrayNull((*boolarray)->blkmem, &(*boolarray)->vals, (*boolarray)->valssize);
   BMSfreeBlockMemory((*boolarray)->blkmem, boolarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPboolarrayExtend(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(boolarray != NULL);
   assert(boolarray->minusedidx == INT_MAX || boolarray->firstidx >= 0);
   assert(boolarray->maxusedidx == INT_MIN || boolarray->firstidx >= 0);
   assert(boolarray->minusedidx == INT_MAX || boolarray->minusedidx >= boolarray->firstidx);
   assert(boolarray->maxusedidx == INT_MIN || boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);
   
   minidx = MIN(minidx, boolarray->minusedidx);
   maxidx = MAX(maxidx, boolarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   SCIPdebugMessage("extending boolarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n", 
      (void*)boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > boolarray->valssize )
   {
      SCIP_Bool* newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = SCIPsetCalcMemGrowSize(set, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(boolarray->blkmem, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( boolarray->firstidx != -1 )
      {
         for( i = 0; i < boolarray->minusedidx - newfirstidx; ++i )
            newvals[i] = FALSE;

         /* check for possible overflow or negative value */
         assert(boolarray->maxusedidx - boolarray->minusedidx + 1 > 0);

         BMScopyMemoryArray(&newvals[boolarray->minusedidx - newfirstidx],
            &boolarray->vals[boolarray->minusedidx - boolarray->firstidx],
            boolarray->maxusedidx - boolarray->minusedidx + 1); /*lint !e866*/
         for( i = boolarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = FALSE;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = FALSE;
      }

      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(boolarray->blkmem, &boolarray->vals, boolarray->valssize);
      boolarray->vals = newvals;
      boolarray->valssize = newvalssize;
      boolarray->firstidx = newfirstidx;
   }
   else if( boolarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = boolarray->valssize - nused;
      assert(nfree >= 0);
      boolarray->firstidx = minidx - nfree/2;
      assert(boolarray->firstidx <= minidx);
      assert(maxidx < boolarray->firstidx + boolarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < boolarray->valssize; ++i )
         assert(boolarray->vals[i] == FALSE);
#endif
   }
   else if( minidx < boolarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = boolarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + boolarray->valssize);
      
      if( boolarray->minusedidx <= boolarray->maxusedidx )
      {
         int shift;

         assert(boolarray->firstidx <= boolarray->minusedidx);
         assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);

         /* shift used part of array to the right */
         shift = boolarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = boolarray->maxusedidx - boolarray->firstidx; i >= boolarray->minusedidx - boolarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < boolarray->valssize);
            boolarray->vals[i + shift] = boolarray->vals[i];
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            boolarray->vals[boolarray->minusedidx - boolarray->firstidx + i] = FALSE;
      }
      boolarray->firstidx = newfirstidx;
   }
   else if( maxidx >= boolarray->firstidx + boolarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = boolarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + boolarray->valssize);
      
      if( boolarray->minusedidx <= boolarray->maxusedidx )
      {
         int shift;

         assert(boolarray->firstidx <= boolarray->minusedidx);
         assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - boolarray->firstidx;
         assert(shift > 0);
         for( i = boolarray->minusedidx - boolarray->firstidx; i <= boolarray->maxusedidx - boolarray->firstidx; ++i )
         {
            assert(0 <= i - shift && i - shift < boolarray->valssize);
            boolarray->vals[i - shift] = boolarray->vals[i];
         }
         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            boolarray->vals[boolarray->maxusedidx - boolarray->firstidx - i] = FALSE;
      }
      boolarray->firstidx = newfirstidx;
   }

   assert(minidx >= boolarray->firstidx);
   assert(maxidx < boolarray->firstidx + boolarray->valssize);

   return SCIP_OKAY;
}

/** clears a dynamic bool array */
SCIP_RETCODE SCIPboolarrayClear(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(boolarray != NULL);

   SCIPdebugMessage("clearing boolarray %p (firstidx=%d, size=%d, range=[%d,%d])\n", 
      (void*)boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx);

   if( boolarray->minusedidx <= boolarray->maxusedidx )
   {
      int i;
   
      assert(boolarray->firstidx <= boolarray->minusedidx);
      assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
      assert(boolarray->firstidx != -1);
      assert(boolarray->valssize > 0);

      /* clear the used part of array */
      for( i = boolarray->minusedidx - boolarray->firstidx; i <= boolarray->maxusedidx - boolarray->firstidx; ++i )
         boolarray->vals[i] = FALSE;

      /* mark the array cleared */
      boolarray->minusedidx = INT_MAX;
      boolarray->maxusedidx = INT_MIN;
   }
   assert(boolarray->minusedidx == INT_MAX);
   assert(boolarray->maxusedidx == INT_MIN);

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
SCIP_Bool SCIPboolarrayGetVal(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(boolarray != NULL);
   assert(idx >= 0);
   
   if( idx < boolarray->minusedidx || idx > boolarray->maxusedidx )
      return FALSE;
   else
   {
      assert(boolarray->vals != NULL);
      assert(idx - boolarray->firstidx >= 0);
      assert(idx - boolarray->firstidx < boolarray->valssize);

      return boolarray->vals[idx - boolarray->firstidx];
   }
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPboolarraySetVal(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   idx,                /**< array index to set value for */
   SCIP_Bool             val                 /**< value to set array index to */
   )
{
   assert(boolarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting boolarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %u\n", 
      (void*)boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx, idx, val);

   if( val != FALSE )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPboolarrayExtend(boolarray, set, idx, idx) );
      assert(idx >= boolarray->firstidx);
      assert(idx < boolarray->firstidx + boolarray->valssize);
      
      /* set the array value of the index */
      boolarray->vals[idx - boolarray->firstidx] = val;

      /* update min/maxusedidx */
      boolarray->minusedidx = MIN(boolarray->minusedidx, idx);
      boolarray->maxusedidx = MAX(boolarray->maxusedidx, idx);
   }
   else if( idx >= boolarray->firstidx && idx < boolarray->firstidx + boolarray->valssize )
   {
      /* set the array value of the index to zero */
      boolarray->vals[idx - boolarray->firstidx] = FALSE;
      
      /* check, if we can tighten the min/maxusedidx */
      if( idx == boolarray->minusedidx )
      {
         assert(boolarray->maxusedidx >= 0);
         assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
         do
         {
            boolarray->minusedidx++;
         }
         while( boolarray->minusedidx <= boolarray->maxusedidx
            && boolarray->vals[boolarray->minusedidx - boolarray->firstidx] == FALSE );
         if( boolarray->minusedidx > boolarray->maxusedidx )
         {
            boolarray->minusedidx = INT_MAX;
            boolarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == boolarray->maxusedidx )
      {
         assert(boolarray->minusedidx >= 0);
         assert(boolarray->minusedidx < boolarray->maxusedidx);
         assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
         do
         {
            boolarray->maxusedidx--;
            assert(boolarray->minusedidx <= boolarray->maxusedidx);
         }
         while( boolarray->vals[boolarray->maxusedidx - boolarray->firstidx] == FALSE );
      }      
   }

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPboolarrayGetMinIdx(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(boolarray != NULL);

   return boolarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPboolarrayGetMaxIdx(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(boolarray != NULL);

   return boolarray->maxusedidx;
}


/** creates a dynamic array of pointer values */
SCIP_RETCODE SCIPptrarrayCreate(
   SCIP_PTRARRAY**       ptrarray,           /**< pointer to store the ptr array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(ptrarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, ptrarray) );
   (*ptrarray)->blkmem = blkmem;
   (*ptrarray)->vals = NULL;
   (*ptrarray)->valssize = 0;
   (*ptrarray)->firstidx = -1;
   (*ptrarray)->minusedidx = INT_MAX;
   (*ptrarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of pointer values */
SCIP_RETCODE SCIPptrarrayCopy(
   SCIP_PTRARRAY**       ptrarray,           /**< pointer to store the copied ptr array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PTRARRAY*        sourceptrarray      /**< dynamic ptr array to copy */
   )
{
   assert(ptrarray != NULL);
   assert(sourceptrarray != NULL);

   SCIP_CALL( SCIPptrarrayCreate(ptrarray, blkmem) );
   if( sourceptrarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*ptrarray)->vals, sourceptrarray->vals, sourceptrarray->valssize) );
   }
   (*ptrarray)->valssize = sourceptrarray->valssize;
   (*ptrarray)->firstidx = sourceptrarray->firstidx;
   (*ptrarray)->minusedidx = sourceptrarray->minusedidx;
   (*ptrarray)->maxusedidx = sourceptrarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of pointer values */
SCIP_RETCODE SCIPptrarrayFree(
   SCIP_PTRARRAY**       ptrarray            /**< pointer to the ptr array */
   )
{
   assert(ptrarray != NULL);
   assert(*ptrarray != NULL);

   BMSfreeBlockMemoryArrayNull((*ptrarray)->blkmem, &(*ptrarray)->vals, (*ptrarray)->valssize);
   BMSfreeBlockMemory((*ptrarray)->blkmem, ptrarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPptrarrayExtend(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(ptrarray != NULL);
   assert(ptrarray->minusedidx == INT_MAX || ptrarray->firstidx >= 0);
   assert(ptrarray->maxusedidx == INT_MIN || ptrarray->firstidx >= 0);
   assert(ptrarray->minusedidx == INT_MAX || ptrarray->minusedidx >= ptrarray->firstidx);
   assert(ptrarray->maxusedidx == INT_MIN || ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);
   
   minidx = MIN(minidx, ptrarray->minusedidx);
   maxidx = MAX(maxidx, ptrarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   SCIPdebugMessage("extending ptrarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n", 
      (void*)ptrarray, ptrarray->firstidx, ptrarray->valssize, ptrarray->minusedidx, ptrarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > ptrarray->valssize )
   {
      void** newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = SCIPsetCalcMemGrowSize(set, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(ptrarray->blkmem, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( ptrarray->firstidx != -1 )
      {
         for( i = 0; i < ptrarray->minusedidx - newfirstidx; ++i )
            newvals[i] = NULL;

         /* check for possible overflow or negative value */
         assert(ptrarray->maxusedidx - ptrarray->minusedidx + 1 > 0);

         BMScopyMemoryArray(&newvals[ptrarray->minusedidx - newfirstidx],
            &(ptrarray->vals[ptrarray->minusedidx - ptrarray->firstidx]),
            ptrarray->maxusedidx - ptrarray->minusedidx + 1); /*lint !e866*/
         for( i = ptrarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = NULL;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = NULL;
      }
      
      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(ptrarray->blkmem, &ptrarray->vals, ptrarray->valssize);
      ptrarray->vals = newvals;
      ptrarray->valssize = newvalssize;
      ptrarray->firstidx = newfirstidx;
   }
   else if( ptrarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = ptrarray->valssize - nused;
      assert(nfree >= 0);
      ptrarray->firstidx = minidx - nfree/2;
      assert(ptrarray->firstidx <= minidx);
      assert(maxidx < ptrarray->firstidx + ptrarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < ptrarray->valssize; ++i )
         assert(ptrarray->vals[i] == NULL);
#endif
   }
   else if( minidx < ptrarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = ptrarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + ptrarray->valssize);
      
      if( ptrarray->minusedidx <= ptrarray->maxusedidx )
      {
         int shift;

         assert(ptrarray->firstidx <= ptrarray->minusedidx);
         assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);

         /* shift used part of array to the right */
         shift = ptrarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = ptrarray->maxusedidx - ptrarray->firstidx; i >= ptrarray->minusedidx - ptrarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < ptrarray->valssize);
            ptrarray->vals[i + shift] = ptrarray->vals[i];
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            ptrarray->vals[ptrarray->minusedidx - ptrarray->firstidx + i] = NULL;
      }
      ptrarray->firstidx = newfirstidx;
   }
   else if( maxidx >= ptrarray->firstidx + ptrarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = ptrarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + ptrarray->valssize);
      
      if( ptrarray->minusedidx <= ptrarray->maxusedidx )
      {
         int shift;

         assert(ptrarray->firstidx <= ptrarray->minusedidx);
         assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - ptrarray->firstidx;
         assert(shift > 0);
         for( i = ptrarray->minusedidx - ptrarray->firstidx; i <= ptrarray->maxusedidx - ptrarray->firstidx; ++i )
         {
            assert(0 <= i - shift && i - shift < ptrarray->valssize);
            ptrarray->vals[i - shift] = ptrarray->vals[i];
         }
         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            ptrarray->vals[ptrarray->maxusedidx - ptrarray->firstidx - i] = NULL;
      }
      ptrarray->firstidx = newfirstidx;
   }

   assert(minidx >= ptrarray->firstidx);
   assert(maxidx < ptrarray->firstidx + ptrarray->valssize);

   return SCIP_OKAY;
}

/** clears a dynamic pointer array */
SCIP_RETCODE SCIPptrarrayClear(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(ptrarray != NULL);

   SCIPdebugMessage("clearing ptrarray %p (firstidx=%d, size=%d, range=[%d,%d])\n", 
      (void*)ptrarray, ptrarray->firstidx, ptrarray->valssize, ptrarray->minusedidx, ptrarray->maxusedidx);

   if( ptrarray->minusedidx <= ptrarray->maxusedidx )
   {
      int i;
   
      assert(ptrarray->firstidx <= ptrarray->minusedidx);
      assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
      assert(ptrarray->firstidx != -1);
      assert(ptrarray->valssize > 0);

      /* clear the used part of array */
      for( i = ptrarray->minusedidx - ptrarray->firstidx; i <= ptrarray->maxusedidx - ptrarray->firstidx; ++i )
         ptrarray->vals[i] = NULL;

      /* mark the array cleared */
      ptrarray->minusedidx = INT_MAX;
      ptrarray->maxusedidx = INT_MIN;
   }
   assert(ptrarray->minusedidx == INT_MAX);
   assert(ptrarray->maxusedidx == INT_MIN);

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
void* SCIPptrarrayGetVal(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   idx                 /**< array index to get value for */
   )
{
   assert(ptrarray != NULL);
   assert(idx >= 0);
   
   if( idx < ptrarray->minusedidx || idx > ptrarray->maxusedidx )
      return NULL;
   else
   {
      assert(ptrarray->vals != NULL);
      assert(idx - ptrarray->firstidx >= 0);
      assert(idx - ptrarray->firstidx < ptrarray->valssize);

      return ptrarray->vals[idx - ptrarray->firstidx];
   }
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPptrarraySetVal(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   idx,                /**< array index to set value for */
   void*                 val                 /**< value to set array index to */
   )
{
   assert(ptrarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting ptrarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %p\n", 
      (void*)ptrarray, ptrarray->firstidx, ptrarray->valssize, ptrarray->minusedidx, ptrarray->maxusedidx, idx, val);

   if( val != NULL )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPptrarrayExtend(ptrarray, set, idx, idx) );
      assert(idx >= ptrarray->firstidx);
      assert(idx < ptrarray->firstidx + ptrarray->valssize);
      
      /* set the array value of the index */
      ptrarray->vals[idx - ptrarray->firstidx] = val;

      /* update min/maxusedidx */
      ptrarray->minusedidx = MIN(ptrarray->minusedidx, idx);
      ptrarray->maxusedidx = MAX(ptrarray->maxusedidx, idx);
   }
   else if( idx >= ptrarray->firstidx && idx < ptrarray->firstidx + ptrarray->valssize )
   {
      /* set the array value of the index to zero */
      ptrarray->vals[idx - ptrarray->firstidx] = NULL;
      
      /* check, if we can tighten the min/maxusedidx */
      if( idx == ptrarray->minusedidx )
      {
         assert(ptrarray->maxusedidx >= 0);
         assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
         do
         {
            ptrarray->minusedidx++;
         }
         while( ptrarray->minusedidx <= ptrarray->maxusedidx
            && ptrarray->vals[ptrarray->minusedidx - ptrarray->firstidx] == NULL );
         if( ptrarray->minusedidx > ptrarray->maxusedidx )
         {
            ptrarray->minusedidx = INT_MAX;
            ptrarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == ptrarray->maxusedidx )
      {
         assert(ptrarray->minusedidx >= 0);
         assert(ptrarray->minusedidx < ptrarray->maxusedidx);
         assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
         do
         {
            ptrarray->maxusedidx--;
            assert(ptrarray->minusedidx <= ptrarray->maxusedidx);
         }
         while( ptrarray->vals[ptrarray->maxusedidx - ptrarray->firstidx] == NULL );
      }      
   }

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPptrarrayGetMinIdx(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(ptrarray != NULL);

   return ptrarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPptrarrayGetMaxIdx(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(ptrarray != NULL);

   return ptrarray->maxusedidx;
}




/*
 * Sorting algorithms
 */

/* first all upwards-sorting methods */

/** sort an indexed element set in non-decreasing order, resulting in a permutation index array */
void SCIPsort(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   )
{
   int pos;

   assert(indcomp != NULL);
   assert(len == 0 || perm != NULL);

   /* create identity permutation */
   for( pos = 0; pos < len; ++pos )
      perm[pos] = pos;
   
   SCIPsortInd(perm, indcomp, dataptr, len);
}

/* SCIPsortInd(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Ind
#define SORTTPL_KEYTYPE     int
#define SORTTPL_INDCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Ptr
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtr
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrReal
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrRealInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrReal
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrRealIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrRealInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrLongInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrLongInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrLongIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrLongIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Real
#define SORTTPL_KEYTYPE     SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealPtrPtrInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealIntLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealIntLong
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealIntPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealLongRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealLongRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Longint
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealIntInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealRealRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealRealRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealRealRealBoolPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealRealBoolPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Int
#define SORTTPL_KEYTYPE     int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntLong
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntPtrReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntPtrIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntPtrIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Long
#define SORTTPL_KEYTYPE     SCIP_Longint
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtr
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrPtrInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrPtrIntInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrPtrBoolInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrPtrBoolInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/



/* now all downwards-sorting methods */


/** sort an indexed element set in non-increasing order, resulting in a permutation index array */
void SCIPsortDown(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   )
{
   int pos;

   assert(indcomp != NULL);
   assert(len == 0 || perm != NULL);

   /* create identity permutation */
   for( pos = 0; pos < len; ++pos )
      perm[pos] = pos;
   
   SCIPsortDownInd(perm, indcomp, dataptr, len);
}


/* SCIPsortDownInd(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownInd
#define SORTTPL_KEYTYPE     int
#define SORTTPL_INDCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtr
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtr
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrReal
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortDownPtrBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortDownPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrRealInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrReal
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrRealIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrLongInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrLongInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrLongIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrLongIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownReal
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealIntLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealIntLong
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealIntPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealPtrPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealPtrPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealLongRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealLongRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Longint
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealIntInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealRealBoolPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealRealBoolPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortDownIntIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortDownIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntLong
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntPtrIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntPtrIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  SCIP_Real
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLong
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtr
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrPtrInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrPtrIntInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrPtrBoolInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrPtrBoolInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/*
 * Stair map
 */

/** creates stair map */
SCIP_RETCODE SCIPstairmapCreate(
   SCIP_STAIRMAP**       stairmap,           /**< pointer to store the created stair map */
   int                   upperbound,         /**< upper bound of the stairmap */
   int                   ntimepoints         /**< minimum size to ensure */
   )
{
   assert(stairmap != NULL);
   assert(upperbound > 0);
   assert(ntimepoints > 0);

   SCIP_ALLOC( BMSallocMemory(stairmap) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*stairmap)->timepoints, ntimepoints) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*stairmap)->freecapacities, ntimepoints) );

   /* setup cumulative stairmap for use */
   (*stairmap)->ntimepoints = 2;
   (*stairmap)->timepoints[0] = 0;
   (*stairmap)->timepoints[1] = INT_MAX;
   (*stairmap)->freecapacities[0] = upperbound;
   (*stairmap)->freecapacities[1] = 0;
   (*stairmap)->arraysize = ntimepoints;
   
   return SCIP_OKAY;
}

/** frees given stair map */
void SCIPstairmapFree(
   SCIP_STAIRMAP**       stairmap            /**< pointer to the stair map */
   )
{
   assert(stairmap != NULL);
   assert(*stairmap != NULL);
   
   /* free main hash map data structure */
   BMSfreeMemoryArray(&(*stairmap)->freecapacities);
   BMSfreeMemoryArray(&(*stairmap)->timepoints);
   BMSfreeMemory(stairmap);
}

/** resizes the stair map arrays */
SCIP_RETCODE SCIPstairmapResize(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to resize */
   int                   ntimepoints         /**< minimum size to ensure */
   )
{
   assert(stairmap != NULL);
   assert(ntimepoints >= 0);
   assert(stairmap->timepoints != NULL);
   assert(stairmap->freecapacities != NULL);
   
   if( stairmap->ntimepoints >= ntimepoints )
      return SCIP_OKAY;

   /* grow arrays of times and free capacity */
   SCIP_ALLOC( BMSreallocMemoryArray(&stairmap->timepoints, ntimepoints) );
   SCIP_ALLOC( BMSreallocMemoryArray(&stairmap->freecapacities, ntimepoints) );
   stairmap->arraysize = ntimepoints;

   return SCIP_OKAY;
}

/** output of the given stair map */
void SCIPstairmapPrint(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to output */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int t;
   
   for( t = 0; t < stairmap->ntimepoints; ++t )
   {
      SCIPmessageFPrintInfo(file, "i: %d, tp: %d, fc: %d ;", t, stairmap->timepoints[t], stairmap-> freecapacities[t]); 
   }
   
   SCIPmessageFPrintInfo(file,"\n");
}

/** returns if the given time point exists in the stair map and stores the position of the given time point if it
 *  exists; otherwise the position of the next smaller existing time point is stored
 */
static
SCIP_Bool stairmapFindLeft(
   SCIP_STAIRMAP*        stairmap,             /**< stair map to search */
   int                   timepoint,            /**< time point to search for */
   int*                  pos                   /**< pointer to store the position */
   )
{
   assert(stairmap != NULL);
   assert(timepoint >= 0);
   assert(stairmap->ntimepoints > 0);
   assert(stairmap->timepoints[0] == 0);

   /* find the position of time point in the time points array via binary search */
   if( SCIPsortedvecFindInt(stairmap->timepoints, timepoint, stairmap->ntimepoints, pos) )
      return TRUE;
 
   assert(*pos > 0);
   (*pos)--;
   
   return FALSE;
}

/** inserts the given time point into the stairmap if it this time point does not exists yet; returns its position in the
 *  time point array 
 */
static
int stairmapInsertTimepoint(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to insert the time point */
   int                   timepoint           /**< time point to insert */
   )
{
   int pos;

   assert(stairmap != NULL);
   assert(timepoint >= 0);
   assert(stairmap->arraysize >= stairmap->ntimepoints);

   if( timepoint == 0 )
      return 0;
      
   /* get the position of the given time point in the stair map array if it exists; otherwise the position of the next
    * smaller existing time point 
    */
   if( stairmapFindLeft(stairmap, timepoint, &pos) )
   {
      /* if the time point exists return the corresponding position */
      assert(pos >= 0 && pos < stairmap->ntimepoints);
      return pos;
   }

   assert(pos >= 0 && pos < stairmap->ntimepoints);
   assert(timepoint >= stairmap->timepoints[pos]);
   assert(pos + 1 < stairmap->arraysize);

   /* insert new time point into the (sorted) stair map */
   SCIPsortedvecInsertIntInt(stairmap->timepoints, stairmap->freecapacities, timepoint, stairmap->freecapacities[pos], 
      &stairmap->ntimepoints);
   
#ifndef NDEBUG
   /* check if the time points are sorted */
   {
      int i;   
      for( i = 1; i < stairmap->ntimepoints; ++i )
         assert(stairmap->timepoints[i-1] < stairmap->timepoints[i]);
   }
#endif
   
   return pos+1;
}

/** updates the stair map due to inserting of a core */
static
void stairmapUpdate(
   SCIP_STAIRMAP*        stairmap,           /**< stairmap to update */
   int                   left,               /**< left side of core interval */
   int                   right,              /**< right side of core interval */
   int                   height,             /**< height of the core */
   SCIP_Bool*            infeasible          /**< pointer to store if the update is infeasible */
   )
{
   int startpos;
   int endpos;
   int i;

   assert(stairmap != NULL);
   assert(stairmap->arraysize >= stairmap->ntimepoints);
   assert(left >= 0);
   assert(left < right);
   assert(infeasible != NULL);
   
   (*infeasible) = FALSE;
   
   /* get position of the starttime in stairmap */
   startpos = stairmapInsertTimepoint(stairmap, left);
   assert(stairmap->timepoints[startpos] == left);

   /* get position of the endtime in stairmap */
   endpos = stairmapInsertTimepoint(stairmap, right);
   assert(stairmap->timepoints[endpos] == right);

   assert(startpos < endpos);
   assert(stairmap->arraysize >= stairmap->ntimepoints);

   /* remove/add the given height from the stair map */
   for( i = startpos; i < endpos; ++i )
   {
      stairmap->freecapacities[i] -= height;

      if( stairmap->freecapacities[i] < 0 )
      {
         *infeasible = TRUE;

         /* remove infeasible core */
         for( ; i >= startpos; --i ) /*lint !e445*/
            stairmap->freecapacities[i] += height;
         
         break;
      }      
   }
}

/** insert a core into stair map; if core is non-empty the stair map will be updated otherwise nothing happens */
void SCIPstairmapInsertCore(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to use */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   height,             /**< height of the core */
   SCIP_Bool*            infeasible          /**< pointer to store if the core does not fit due to capacity */
   )
{
   assert(stairmap != NULL);
   assert(left < right);
   assert(infeasible != NULL);
   
   (*infeasible) = FALSE;
   
   /* insert core into the stair map */
   SCIPdebugMessage("insert core [%d,%d] with height %d\n", left, right, height);
   
   /* try to insert core into the stair map */
   stairmapUpdate(stairmap, left, right, height, infeasible);
}

/** subtracts the height from the stair map during core time */
void SCIPstairmapDeleteCore(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to use */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   height              /**< height of the core */
   )
{
   SCIP_Bool infeasible;
   
   assert(left < right);
#ifndef NDEBUG
      {
         /* check if the begin and end time points of the core correspond to a time point in the stairmap; this should be
          * the case since we added the core before to the stair map 
          */
         int pos;
         assert(stairmapFindLeft(stairmap, left, &pos));
         assert(stairmapFindLeft(stairmap, right, &pos));
      }
#endif
      
      /* remove the core from the current stair map */
      SCIPdebugMessage("delete core [%d,%d] with height %d\n", left, right, height);
      
      stairmapUpdate(stairmap, left, right, -height, &infeasible);
      assert(!infeasible);
}
   
/** returns the time point at the given position */   
int SCIPstairmapGetTimepoint(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to use */
   int                   pos                 /**< position */
   )
{
   assert(stairmap != NULL);
   assert(pos < stairmap->ntimepoints);

   return stairmap->timepoints[pos];
}

/** returns TRUE if the core  (given by its height and during) can be inserted at the given time point; otherwise FALSE */
SCIP_Bool SCIPstairmapIsFeasibleStart(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to use */
   int                   timepoint,          /**< time point to start */
   int                   duration,           /**< duration of the core */
   int                   height,             /**< height of the core */
   int*                  pos                 /**< pointer to store the earliest position where the core does not fit */
   )
{
   int endtime;
   int startpos;
   int endpos;
   int p;

   assert(stairmap != NULL);
   assert(timepoint >= 0);
   assert(height >= 0);
   assert(pos != NULL);

   if( duration == 0 )
      return TRUE;
   
   endtime = timepoint + duration; 

   /* check if the activity fits at timepoint */
   (void)stairmapFindLeft(stairmap, timepoint, &startpos);

   if( !stairmapFindLeft(stairmap, endtime, &endpos) )
      endpos++;
   
   assert(stairmap->timepoints[startpos] <= timepoint);
   assert(stairmap->timepoints[endpos] >= endtime);
   
   for( p = startpos; p < endpos; ++p )
   {
      if( stairmap->freecapacities[p] < height )
      {
         (*pos) = p;
         return FALSE;
      }
   }
   
   return TRUE;
}

/** return the earliest possible starting point within the time interval [lb,ub] for a given core (given by its height
 *  and duration)
 */
int SCIPstairmapGetEarliestFeasibleStart(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to use */
   int                   lb,                 /**< earliest starting time of the given core */
   int                   ub,                 /**< latest starting time of the given core */
   int                   duration,           /**< duration of the core */
   int                   height,             /**< height of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the core cannot be inserted */
   )
{
   int starttime;
   int pos;
   
   assert(stairmap != NULL);
   assert(lb >= 0);
   assert(duration >= 0);
   assert(height >= 0);
   assert(infeasible != NULL);
   assert(stairmap->timepoints[stairmap->ntimepoints-1] > ub);

   if( lb > ub )
   {
      *infeasible = TRUE;
      return lb;
   }
   
   if( duration == 0 || height == 0 ) 
   {
      *infeasible = FALSE;
      return lb;
   }

   starttime = lb;

   (void)stairmapFindLeft(stairmap, starttime, &pos);
   assert(stairmap->timepoints[pos] <= starttime);
   
   (*infeasible) = TRUE;

   while( (*infeasible) && starttime <= ub )
   {
      if( SCIPstairmapIsFeasibleStart(stairmap, starttime, duration, height, &pos) )
      {
         (*infeasible) = FALSE;
         return starttime;
      }
    
      /* the core did not fit into the stair map since at time point "pos" not enough capacity is available; therefore we
       * can proceed with the next time point
       */
      assert(stairmap->freecapacities[pos] < height);
      pos++;
      
      /* check if we exceed the time point array */
      if( pos >= stairmap->ntimepoints )
         break;
      
      starttime = stairmap->timepoints[pos];
   }
   
   assert(*infeasible || starttime <= ub);
   return starttime;
}

/** return the latest possible starting point within the time interval [lb,ub] for a given core (given by its height and
 *  duration)
 */
int SCIPstairmapGetLatestFeasibleStart(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to use */
   int                   lb,                 /**< earliest possible start point */
   int                   ub,                 /**< latest possible start point */
   int                   duration,           /**< duration of the core */
   int                   height,             /**< height of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the core cannot be inserted */
   )
{
   int starttime;
   int pos;

   assert(stairmap != NULL);
   assert(lb >= 0);
   assert(lb <= ub);
   assert(duration >= 0);
   assert(height >= 0);
   assert(infeasible != NULL);
   assert(stairmap->timepoints[stairmap->ntimepoints-1] > ub);
   
   if( duration == 0 || height == 0 ) 
      return ub;

   starttime = ub;   
   (void)stairmapFindLeft(stairmap, starttime, &pos);
   assert(stairmap->timepoints[pos] <= starttime);
   
   (*infeasible) = TRUE;

   while( (*infeasible) && starttime >= lb )
   {
      if( SCIPstairmapIsFeasibleStart(stairmap, starttime, duration, height, &pos) )
      {
         (*infeasible) = FALSE;
         return starttime;
      }
      assert(pos >= 0);

      /* the core did not fit into the stair map since at time point "pos" not enough capacity is available; 
       * therefore we can proceed with the next time point  */
      assert(stairmap->freecapacities[pos] < height);
            
      starttime = stairmap->timepoints[pos] - duration;
   }

   assert(*infeasible || starttime >= lb);
  
   return starttime;
}

/*
 * Numerical methods
 */

/** returns the machine epsilon: the smallest number eps > 0, for which 1.0 + eps > 1.0 */
SCIP_Real SCIPcalcMachineEpsilon(
   void
   )
{
   SCIP_Real eps;
   SCIP_Real lasteps;
   SCIP_Real one;
   SCIP_Real onepluseps;

   one = 1.0;
   eps = 1.0;
   do
   {
      lasteps = eps;
      eps /= 2.0;
      onepluseps = one + eps;
   }
   while( onepluseps > one );

   return lasteps;
}

/** calculates the greatest common divisor of the two given values */
SCIP_Longint SCIPcalcGreComDiv(
   SCIP_Longint          val1,               /**< first value of greatest common devisor calculation */
   SCIP_Longint          val2                /**< second value of greatest common devisor calculation */
   )
{
   int t;

   assert(val1 > 0);
   assert(val2 > 0);
   
   t = 0;
   /* if val1 is even, divide it by 2 */
   while( !(val1 & 1) )
   {
      val1 >>= 1; /*lint !e704*/
      
      /* if val2 is even too, divide it by 2 and increase t(=number of e) */
      if( !(val2 & 1) )
      {
         val2 >>= 1; /*lint !e704*/
         ++t;
      }
      /* only val1 can be odd */
      else
      {
         /* while val1 is even, divide it by 2 */
         while( !(val1 & 1) )
            val1 >>= 1; /*lint !e704*/
         
         break;
      }
   }
   /* while val2 is even, divide it by 2 */
   while( !(val2 & 1) )
      val2 >>= 1; /*lint !e704*/
   
   /* val1 and val 2 are odd */
   while( val1 != val2 )
   {
      if( val1 > val2 )
      {
         val1 -= val2;
         /* val1 is now even, divide it by 2  */
         do 
         {
            val1 >>= 1;   /*lint !e704*/
         } while( !(val1 & 1) );
      }
      else 
      {
         val2 -= val1;
         /* val2 is now even, divide it by 2  */
         do 
         {
            val2 >>= 1;  /*lint !e704*/
         } while( !(val2 & 1) );
      }
   }

   return (val1 << t);  /*lint !e703*/
}

/** calculates the smallest common multiple of the two given values */
SCIP_Longint SCIPcalcSmaComMul(
   SCIP_Longint          val1,               /**< first value of smallest common multiple calculation */
   SCIP_Longint          val2                /**< second value of smallest common multiple calculation */
   )
{
   SCIP_Longint gcd;

   assert(val1 > 0);
   assert(val2 > 0);

   gcd = SCIPcalcGreComDiv(val1, val2);
   
   return val1/gcd * val2;
}

static const SCIP_Real simplednoms[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
                                   17.0, 18.0, 19.0, 25.0, -1.0};

/** converts a real number into a (approximate) rational representation, and returns TRUE iff the conversion was
 *  successful
 */
SCIP_Bool SCIPrealToRational(
   SCIP_Real             val,                /**< real value r to convert into rational number */
   SCIP_Real             mindelta,           /**< minimal allowed difference r - q of real r and rational q = n/d */
   SCIP_Real             maxdelta,           /**< maximal allowed difference r - q of real r and rational q = n/d */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   )
{
   SCIP_Real a;
   SCIP_Real b;
   SCIP_Real g0;
   SCIP_Real g1;
   SCIP_Real gx;
   SCIP_Real h0;
   SCIP_Real h1;
   SCIP_Real hx;
   SCIP_Real delta0;
   SCIP_Real delta1;
   SCIP_Real epsilon;
   int i;

   assert(mindelta < 0.0);
   assert(maxdelta > 0.0);
   assert(nominator != NULL);
   assert(denominator != NULL);

   /* try the simple denominators first: each value of the simpledenoms table multiplied by powers of 10
    * is tried as denominator
    */
   for( i = 0; simplednoms[i] > 0.0; ++i )
   {
      SCIP_Real nom;
      SCIP_Real dnom;
      SCIP_Real ratval0;
      SCIP_Real ratval1;

      /* try powers of 10 (including 10^0) */
      dnom = simplednoms[i];
      while( dnom <= maxdnom )
      {
         nom = floor(val * dnom);
         ratval0 = nom/dnom;
         ratval1 = (nom+1.0)/dnom;
         if( mindelta <= val - ratval0 && val - ratval1 <= maxdelta )
         {
            if( val - ratval0 <= maxdelta )
            {
               *nominator = (SCIP_Longint)nom;
               *denominator = (SCIP_Longint)dnom;
               return TRUE;
            }
            if( mindelta <= val - ratval1 )
            {
               *nominator = (SCIP_Longint)(nom+1.0);
               *denominator = (SCIP_Longint)dnom;
               return TRUE;
            }
         }
         dnom *= 10.0;
      }
   }

   /* the simple denominators didn't work: calculate rational representation with arbitrary denominator */
   epsilon = MIN(-mindelta, maxdelta)/2.0;

   b = val;
   a = EPSFLOOR(b, epsilon);
   g0 = a;
   h0 = 1.0;
   g1 = 1.0;
   h1 = 0.0;
   delta0 = val - g0/h0;
   delta1 = (delta0 < 0.0 ? val - (g0-1.0)/h0 : val - (g0+1.0)/h0);
  
   while( (delta0 < mindelta || delta0 > maxdelta) && (delta1 < mindelta || delta1 > maxdelta) )
   {
      assert(EPSGT(b, a, epsilon));
      assert(h0 >= 0.0);
      assert(h1 >= 0.0);

      b = 1.0 / (b - a);
      a = EPSFLOOR(b, epsilon);

      assert(a >= 0.0);
      gx = g0;
      hx = h0;

      g0 = a * g0 + g1;
      h0 = a * h0 + h1;

      g1 = gx;
      h1 = hx;
      
      if( h0 > maxdnom )
         return FALSE;
      
      delta0 = val - g0/h0;
      delta1 = (delta0 < 0.0 ? val - (g0-1.0)/h0 : val - (g0+1.0)/h0);
   }

   if( REALABS(g0) > (SCIP_LONGINT_MAX >> 4) || h0 > (SCIP_LONGINT_MAX >> 4) )
      return FALSE;

   assert(h0 > 0.5);

   if( delta0 < mindelta )
   {
      assert(mindelta <= delta1 && delta1 <= maxdelta);
      *nominator = (SCIP_Longint)(g0 - 1.0);
      *denominator = (SCIP_Longint)h0;
   }
   else if( delta0 > maxdelta )
   {
      assert(mindelta <= delta1 && delta1 <= maxdelta);
      *nominator = (SCIP_Longint)(g0 + 1.0);
      *denominator = (SCIP_Longint)h0;
   }
   else
   {
      *nominator = (SCIP_Longint)g0;
      *denominator = (SCIP_Longint)h0;
   }
   assert(*denominator >= 1);
   assert(val - (SCIP_Real)(*nominator)/(SCIP_Real)(*denominator) >= mindelta);
   assert(val - (SCIP_Real)(*nominator)/(SCIP_Real)(*denominator) <= maxdelta);

   return TRUE;
}

/** checks, whether the given scalar scales the given value to an integral number with error in the given bounds */
static
SCIP_Bool isIntegralScalar(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real downval;
   SCIP_Real upval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = floor(sval);
   upval = ceil(sval);

   return (SCIPrelDiff(sval, downval) <= maxdelta || SCIPrelDiff(sval, upval) >= mindelta);
}

/** additional scalars that are tried in integrality scaling */
static const SCIP_Real scalars[] = {3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0};
static const int nscalars = 9;

/** tries to find a value, such that all given values, if scaled with this value become integral */
SCIP_RETCODE SCIPcalcIntegralScalar(
   SCIP_Real*            vals,               /**< values to scale */
   int                   nvals,              /**< number of values to scale */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal allowed scalar */
   SCIP_Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   SCIP_Bool*            success             /**< stores whether returned value is valid */
   )
{
   SCIP_Longint gcd;
   SCIP_Longint scm;
   SCIP_Longint nominator;
   SCIP_Longint denominator;
   SCIP_Real val;
   SCIP_Real minval;
   SCIP_Real absval;
   SCIP_Real scaleval;
   SCIP_Real twomultval;
   SCIP_Bool scalable;
   SCIP_Bool twomult;
   SCIP_Bool rational;
   int c;
   int s;

   assert(vals != NULL);
   assert(nvals >= 0);
   assert(maxdnom >= 1);
   assert(mindelta < 0.0);
   assert(maxdelta > 0.0);
   assert(success != NULL);

   SCIPdebugMessage("trying to find rational representation for given values\n");

   if( intscalar != NULL )
      *intscalar = SCIP_INVALID;
   *success = FALSE;

   /* get minimal absolute non-zero value */
   minval = SCIP_REAL_MAX;
   for( c = 0; c < nvals; ++c )
   {
      val = vals[c];
      if( val < mindelta || val > maxdelta )
      {
         absval = REALABS(val);
         minval = MIN(minval, absval);
      }
   }

   if( minval == SCIP_REAL_MAX )
   {
      /* all coefficients are zero (inside tolerances) */
      if( intscalar != NULL )
         *intscalar = 1.0;
      *success = TRUE;
      SCIPdebugMessage(" -> all values are zero (inside tolerances)\n");

      return SCIP_OKAY;
   }
   assert(minval > MIN(-mindelta, maxdelta));

   /* try, if values can be made integral multiplying them with the reciprocal of the smallest value and a power of 2 */
   scalable = TRUE;
   scaleval = 1.0/minval;
   for( c = 0; c < nvals && scalable; ++c )
   {
      /* check, if the value can be scaled with a simple scalar */
      val = vals[c];
      if( val == 0.0 ) /* zeros are allowed in the vals array */
         continue;
    
      absval = REALABS(val);
      while( scaleval <= maxscale
         && (absval * scaleval < 0.5 || !isIntegralScalar(val, scaleval, mindelta, maxdelta)) )
      {
         for( s = 0; s < nscalars; ++s )
         {
            if( isIntegralScalar(val, scaleval * scalars[s], mindelta, maxdelta) )
            {
               scaleval *= scalars[s];
               break;
            }
         }
         if( s >= nscalars )
            scaleval *= 2.0;
      }
      scalable = (scaleval <= maxscale);
      SCIPdebugMessage(" -> val=%g, scaleval=%g, val*scaleval=%g, scalable=%u\n", 
         val, scaleval, val*scaleval, scalable);
   }
   if( scalable )
   {
      /* make values integral by dividing them by the smallest value (and multiplying them with a power of 2) */
      assert(scaleval <= maxscale);
      if( intscalar != NULL )
         *intscalar = scaleval;
      *success = TRUE;
      SCIPdebugMessage(" -> integrality can be achieved by scaling with %g (minval=%g)\n", scaleval, minval);
      
      return SCIP_OKAY;
   }

   /* try, if values can be made integral by multiplying them by a power of 2 */
   twomult = TRUE;
   twomultval = 1.0;
   for( c = 0; c < nvals && twomult; ++c )
   {
      /* check, if the value can be scaled with a simple scalar */
      val = vals[c];
      if( val == 0.0 ) /* zeros are allowed in the vals array */
         continue;
      
      absval = REALABS(val);
      while( twomultval <= maxscale
         && (absval * twomultval < 0.5 || !isIntegralScalar(val, twomultval, mindelta, maxdelta)) )
      {
         for( s = 0; s < nscalars; ++s )
         {
            if( isIntegralScalar(val, twomultval * scalars[s], mindelta, maxdelta) )
            {
               twomultval *= scalars[s];
               break;
            }
         }
         if( s >= nscalars )
            twomultval *= 2.0;
      }
      twomult = (twomultval <= maxscale);
      SCIPdebugMessage(" -> val=%g, twomult=%g, val*twomult=%g, twomultable=%u\n",
         val, twomultval, val*twomultval, twomult);
   }
   if( twomult )
   {
      /* make values integral by multiplying them with a power of 2 */
      assert(twomultval <= maxscale);
      if( intscalar != NULL )
         *intscalar = twomultval;
      *success = TRUE;
      SCIPdebugMessage(" -> integrality can be achieved by scaling with %g (power of 2)\n", twomultval);
      
      return SCIP_OKAY;
   }

   /* convert each value into a rational number, calculate the greatest common divisor of the nominators
    * and the smallest common multiple of the denominators
    */
   gcd = 1;
   scm = 1;
   rational = TRUE;

   /* first value (to initialize gcd) */
   for( c = 0; c < nvals && rational; ++c )
   {
      val = vals[c];
      if( val == 0.0 ) /* zeros are allowed in the vals array */
         continue;

      rational = SCIPrealToRational(val, mindelta, maxdelta, maxdnom, &nominator, &denominator);
      if( rational && nominator != 0 )
      {
         assert(denominator > 0);
         gcd = ABS(nominator);
         scm = denominator;
         rational = ((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
         SCIPdebugMessage(" -> c=%d first rational: val: %g == %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", gcd=%"SCIP_LONGINT_FORMAT", scm=%"SCIP_LONGINT_FORMAT", rational=%u\n",
            c, val, nominator, denominator, gcd, scm, rational);
         break;
      }
   }

   /* remaining values */
   for( ++c; c < nvals && rational; ++c )
   {
      val = vals[c];
      if( val == 0.0 ) /* zeros are allowed in the vals array */
         continue;

      rational = SCIPrealToRational(val, mindelta, maxdelta, maxdnom, &nominator, &denominator);
      if( rational && nominator != 0 )
      {
         assert(denominator > 0);
         gcd = SCIPcalcGreComDiv(gcd, ABS(nominator));
         scm *= denominator / SCIPcalcGreComDiv(scm, denominator);
         rational = ((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
         SCIPdebugMessage(" -> c=%d next rational : val: %g == %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", gcd=%"SCIP_LONGINT_FORMAT", scm=%"SCIP_LONGINT_FORMAT", rational=%u\n",
            c, val, nominator, denominator, gcd, scm, rational);
      }
      else
      {
         SCIPdebugMessage(" -> failed to convert %g into a rational representation\n", val);
      }
   }

   if( rational )
   {
      /* make values integral by multiplying them with the smallest common multiple of the denominators */
      assert((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
      if( intscalar != NULL )
         *intscalar = (SCIP_Real)scm/(SCIP_Real)gcd;
      *success = TRUE;
      SCIPdebugMessage(" -> integrality can be achieved by scaling with %g (rational:%"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT")\n", 
         (SCIP_Real)scm/(SCIP_Real)gcd, scm, gcd);
   }

   return SCIP_OKAY;
}

/** given a (usually very small) interval, tries to find a rational number with simple denominator (i.e. a small
 *  number, probably multiplied with powers of 10) out of this interval; returns TRUE iff a valid rational
 *  number inside the interval was found
 */
SCIP_Bool SCIPfindSimpleRational(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed for resulting rational number */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   )
{
   SCIP_Real center;
   SCIP_Real delta;

   assert(lb <= ub);

   center = 0.5*(lb+ub);

   /* in order to compute a rational number that is exactly within the bounds (as the user expects),
    * we computed the allowed delta with downward rounding, if available
    */
   if( SCIPintervalHasRoundingControl() )
   {
      SCIP_ROUNDMODE roundmode;

      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      delta = 0.5*(ub-lb);

      SCIPintervalSetRoundingMode(roundmode);
   }
   else
   {
      delta = 0.5*(ub-lb);
   }

   return SCIPrealToRational(center, -delta, +delta, maxdnom, nominator, denominator);
}

/** given a (usually very small) interval, selects a value inside this interval; it is tried to select a rational number
 *  with simple denominator (i.e. a small number, probably multiplied with powers of 10);
 *  if no valid rational number inside the interval was found, selects the central value of the interval
 */
SCIP_Real SCIPselectSimpleValue(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom             /**< maximal denominator allowed for resulting rational number */
   )
{
   SCIP_Real val;

   val = 0.5*(lb+ub);
   if( lb < ub )
   {
      SCIP_Longint nominator;
      SCIP_Longint denominator;
      SCIP_Bool success;
      
      /* try to find a "simple" rational number inside the interval */
      SCIPdebugMessage("simple rational in [%.9f,%.9f]:", lb, ub);
      success = SCIPfindSimpleRational(lb, ub, maxdnom, &nominator, &denominator);
      if( success )
      {
         val = (SCIP_Real)nominator/(SCIP_Real)denominator;
         SCIPdebugPrintf(" %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT" == %.9f\n", nominator, denominator, val);

         if( val - lb < 0.0 || val - ub > 0.0 )
         {
            SCIPdebugPrintf(" value is out of interval bounds by %g -> failed\n", MAX(lb-val, val-ub));
            val = 0.5*(lb+ub);
         }
      }
      else
      {
         SCIPdebugPrintf(" failed\n");
      }
   }
   
   return val;
}




/*
 * Random Numbers
 */

#ifdef NO_RAND_R

#define SCIP_RAND_MAX 32767
/** returns a random number between 0 and SCIP_RAND_MAX */
static
int getRand(
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   SCIP_Longint nextseed;

   assert(seedp != NULL);

   nextseed = (*seedp) * 1103515245 + 12345;
   *seedp = (unsigned int)nextseed;

   return (int)((unsigned int)(nextseed/(2*(SCIP_RAND_MAX+1))) % (SCIP_RAND_MAX+1));
}

#else

#define SCIP_RAND_MAX RAND_MAX

/** returns a random number between 0 and SCIP_RAND_MAX */
static
int getRand(
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   return rand_r(seedp);
}

#endif

/** returns a random integer between minrandval and maxrandval */
int SCIPgetRandomInt(
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   return minrandval + (int) ((maxrandval - minrandval + 1)*(SCIP_Real)getRand(seedp)/(SCIP_RAND_MAX+1.0));
}

/** returns a random real between minrandval and maxrandval */
SCIP_Real SCIPgetRandomReal(
   SCIP_Real             minrandval,         /**< minimal value to return */
   SCIP_Real             maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   return minrandval + (maxrandval - minrandval)*(SCIP_Real)getRand(seedp)/(SCIP_Real)SCIP_RAND_MAX;
}


/*
 * Permutations / Shuffling
 */

/** swaps the addresses of two pointers */
void SCIPswapPointers(
   void**                pointer1,           /**< first pointer */
   void**                pointer2            /**< second pointer */
   )
{
   void* tmp;

   tmp = *pointer1;
   *pointer1 = *pointer2;
   *pointer2 = tmp;
}

/** randomly shuffles parts of an array using the Fisher-Yates algorithm */
void SCIPpermuteArray(
   void**                array,              /**< array to be shuffled */
   int                   begin,              /**< first index that should be subject to shuffling (0 for whole array) */
   int                   end,                /**< last index that should be subject to shuffling (array size for whole array) */
   unsigned int*         randseed            /**< seed value for the random generator */
   ) 
{
   /* loop backwards through all elements and always swap the current last element to a random position */
   while( end > begin+1 ) 
   {
      int i;
      void* tmp;
      
      end--;
      
      /* get a random position into which the last entry should be shuffled */
      i = SCIPgetRandomInt(begin, end, randseed);

      /* swap the last element and the random element */
      tmp = array[i];
      array[i] = array[end];
      array[end] = tmp;
   }
}

/** draws a random subset of disjoint elements from a given set of disjoint elements;
 *  this implementation is suited for the case that nsubelems is considerably smaller then nelems
 */
SCIP_RETCODE SCIPgetRandomSubset(
   void**                set,                /**< original set, from which elements should be drawn */
   int                   nelems,             /**< number of elements in original set */
   void**                subset,             /**< subset in which drawn elements should be stored */
   int                   nsubelems,          /**< number of elements that should be drawn and stored */
   unsigned int          randseed            /**< seed value for random generator */
   )
{
   int i;
   int j;

   /* if both sets are of equal size, we just copy the array */
   if( nelems == nsubelems)
   {
      BMScopyMemoryArray(subset,set,nelems);
      return SCIP_OKAY;
   }

   /* abort, if size of subset is too big */
   if( nsubelems > nelems )
   {
      SCIPerrorMessage("Cannot create %d-elementary subset of %d-elementary set.\n", nsubelems, nelems);      
      return SCIP_INVALIDDATA;
   }
#ifndef NDEBUG
   for( i = 0; i < nsubelems; i++ ) 
      for( j = 0; j < i; j++ ) 
         assert(set[i] != set[j]);
#endif   

   /* draw each element individually */
   i = 0;
   while( i < nsubelems )
   {
      int r;

      r = SCIPgetRandomInt(0, nelems-1, &randseed);
      subset[i] = set[r];

      /* if we get an element that we already had, we will draw again */
      for( j = 0; j < i; j++ ) 
      {
         if( subset[i] == subset[j] ) 
         {
            --i;
            break;
         }
      }
      ++i;
   }
   return SCIP_OKAY;
}



/*
 * Strings
 */

/** prints an error message containing of the given string followed by a string describing the current system error;
 *  prefers to use the strerror_r method, which is threadsafe; on systems where this method does not exist,
 *  NO_STRERROR_R should be defined (see INSTALL), in this case, srerror is used which is not guaranteed to be
 *  threadsafe (on SUN-systems, it actually is) 
 */
void SCIPprintSysError(
   const char*           message             /**< first part of the error message, e.g. the filename */
   )
{
#ifdef NO_STRERROR_R
   char* buf;
   buf = strerror(errno);
#else
   char buf[SCIP_MAXSTRLEN];
   (void) strerror_r(errno, buf, SCIP_MAXSTRLEN);
   buf[SCIP_MAXSTRLEN - 1] = '\0';
#endif
   SCIPmessagePrintError("%s: %s\n", message, buf);
}

/** extracts tokens from strings - wrapper method for strtok_r() */
char* SCIPstrtok(
   char*                 s,                  /**< string to parse */
   const char*           delim,              /**< delimiters for parsing */
   char**                ptrptr              /**< pointer to working char pointer - must stay the same while parsing */
   )
{
#ifdef NO_STRTOK_R
   return strtok(s, delim);
#else
   return strtok_r(s, delim, ptrptr);
#endif
}

/** translates the given string into a string where symbols ", ', and spaces are escaped with a \ prefix */
void SCIPescapeString(
   char*                 t,                  /**< target buffer to store escaped string */
   int                   bufsize,            /**< size of buffer t */
   const char*           s                   /**< string to transform into escaped string */
   )
{
   int len;
   int i;
   int p;

   assert(t != NULL);
   assert(bufsize > 0);

   len = (int)strlen(s);
   for( p = 0, i = 0; i <= len && p < bufsize; ++i, ++p )
   {
      if( s[i] == ' ' || s[i] == '"' || s[i] == '\'' )
      {
         t[p] = '\\';
         p++;
      }
      if( p < bufsize )
         t[p] = s[i];
   }
   t[bufsize-1] = '\0';
}

/* safe version of snprintf */
int SCIPsnprintf(
   char*                 t,                  /**< target string */
   int                   len,                /**< length of the string to copy */
   const char*           s,                  /**< source string */
   ...                                       /**< further parameters */
   )
{
   va_list ap;
   int n;
   
   assert(t != NULL);
   assert(len > 0);

   va_start(ap, s); /*lint !e826*/
   n = vsnprintf(t, (size_t) len, s, ap);
   va_end(ap);
   if( n < 0 || n >= len )
   {
#ifndef NDEBUG
      if( n < 0 )
         SCIPmessagePrintWarning("vsnprintf returned %d\n",n);
#endif
      t[len-1] = '\0';
      n = len-1;
   }
   return n;
}

/** extract the next token as a double value if it is one; in case no value is parsed the endptr is set to str */
SCIP_Bool SCIPstrToRealValue(
   const char*           str,                /**< string to search */
   SCIP_Real*            value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed */
   )
{
   assert(str != NULL);
   assert(value != NULL);
   assert(endptr != NULL);

   *value = strtod(str, endptr);

   if( *endptr != str && *endptr != NULL )
   {
      SCIPdebugMessage("parsed real value <%g>\n", *value);
      return TRUE;
   }
   *endptr = (char*)str;

   SCIPdebugMessage("failed parseing real value <%s>\n", str);

   return FALSE;
}

/** copies the first size characters between a start and end character of str into token, if no error occured endptr
 *  will point to the position after the read part, otherwise it will point to NULL
 */
void SCIPstrCopySection(
   const char*           str,                /**< string to search */
   char                  startchar,          /**< character which defines the beginning */
   char                  endchar,            /**< character which defines the ending */
   char*                 token,              /**< string to store the copy */
   int                   size,               /**< size of the token char array */
   char**                endptr              /**< pointer to store the final string position if successfully parsed,
                                              *   otherwise str */
   )
{
   const char* copystr;
   int nchars;

   assert(str != NULL);
   assert(token != NULL);
   assert(size > 0);
   assert(endptr != NULL);

   nchars = 0;

   copystr = str;

   /* find starting character */
   while( *str != '\0' && *str != startchar )
      ++str;

   /* did not find start character */
   if( *str == '\0' )
   {
      *endptr = (char*)copystr;
      return;
   }

   /* skip start character */
   ++str;

   /* copy string */
   while( *str != '\0' && *str != endchar && nchars < size-1 )
   {
      assert(nchars < SCIP_MAXSTRLEN);
      token[nchars] = *str;
      nchars++;
      ++str;
   }

   /* add end to token */
   token[nchars] = '\0';

   /* if section was longer than size, we want to reach the end of the parsing section anyway */
   if( nchars == (size-1) )
      while( *str != '\0' && *str != endchar )
         ++str;

   /* did not find end character */
   if( *str == '\0' )
   {
      *endptr = (char*)copystr;
      return;
   }

   /* skip end character */
   ++str;

   SCIPdebugMessage("parsed section <%s>\n", token);

   *endptr = (char*) str;
}

/*
 * File methods
 */

/** returns, whether the given file exists */
SCIP_Bool SCIPfileExists(
   const char*           filename            /**< file name */
   )
{
   FILE* f;

   f = fopen(filename, "r");
   if( f == NULL )
      return FALSE;

   fclose(f);

   return TRUE;
}

/** splits filename into path, name, and extension */
void SCIPsplitFilename(
   char*                 filename,           /**< filename to split; is destroyed (but not freed) during process */
   char**                path,               /**< pointer to store path, or NULL if not needed */
   char**                name,               /**< pointer to store name, or NULL if not needed */
   char**                extension,          /**< pointer to store extension, or NULL if not needed */
   char**                compression         /**< pointer to store compression extension, or NULL if not needed */
   )
{
   char* lastslash;
   char* lastbackslash;
   char* lastdot;

   assert(filename != NULL);

   if( path != NULL )
      *path = NULL;
   if( name != NULL )
      *name = NULL;
   if( extension != NULL )
      *extension = NULL;
   if( compression != NULL )
      *compression = NULL;

   /* treat both slashes '/' and '\' as directory delimiters */
   lastslash = strrchr(filename, '/');
   lastbackslash = strrchr(filename, '\\');
   lastslash = MAX(lastslash, lastbackslash); /*lint !e613*/
   lastdot = strrchr(filename, '.');
   if( lastslash != NULL && lastdot != NULL && lastdot < lastslash ) /* is the last dot belonging to the path? */
      lastdot = NULL;

   /* detect known compression extensions */
#ifdef WITH_ZLIB
   if( lastdot != NULL )
   {
      char* compext;

      compext = lastdot+1;
      if( strcmp(compext, "gz") == 0
        || strcmp(compext, "z") == 0
        || strcmp(compext, "Z") == 0 )
      {
         if( compression != NULL )
            *compression = compext;
         *lastdot = '\0';
      }

      /* find again the last dot in the filename without compression extension */
      lastdot = strrchr(filename, '.');
      if( lastslash != NULL && lastdot != NULL && lastdot < lastslash ) /* is the last dot belonging to the path? */
         lastdot = NULL;
   }
#endif

   if( lastslash == NULL )
   {
      if( name != NULL )
         *name = filename;
   }
   else
   {
      if( path != NULL )
         *path = filename;
      if( name != NULL )
         *name = lastslash+1;
      *lastslash = '\0';
   }

   if( lastdot != NULL )
   {
      if( extension != NULL )
         *extension = lastdot+1;
      *lastdot = '\0';
   }
}




/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPrelDiff

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
SCIP_Real SCIPrelDiff(
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real absval1;
   SCIP_Real absval2;
   SCIP_Real quot;

   absval1 = REALABS(val1);
   absval2 = REALABS(val2);
   quot = MAX3(1.0, absval1, absval2);
   
   return (val1-val2)/quot;
}
