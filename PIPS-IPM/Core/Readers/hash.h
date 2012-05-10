/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* hash table definitions
 *
 * PCx beta-2.0  10/31/96.
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Steve Wright.
 *
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#ifndef HashFile
#define HashFile

typedef struct qpnode *qpListPtr;

typedef struct qpnode {
  int      index;
  char     *entry;
  qpListPtr next;
} qpList;

typedef struct {
  qpListPtr *list;
  int      size;
} qpHashTable;

#ifdef __cplusplus
extern "C" {
#endif
  qpHashTable  *NewHashTable(int size);
  int Insert (qpHashTable * table, char * name, int index );
  int GetIndex( qpHashTable * table, char name[] );
  int DeleteHashTable(qpHashTable * table);
#ifdef __cplusplus
}
#endif

#endif
