/* This file is part of OOPS.
 *
 * OOPS is (c) 2003-2009 Jacek Gondzio and Andreas Grothey, 
 *                       University of Edinburgh
 *
 * OOPS is distributed in a restricted form in the hope that it will be a useful
 * example of what can be done with SML, however it is NOT released under a free
 * software license.
 *
 * You may only redistribute this version of OOPS with a version of SML. You
 * may not link OOPS with code which is not part of SML licensed under the
 * LGPL v3.
 *
 * You may NOT modify, disassemble, or otherwise reverse engineer OOPS.
 *
 * OOPS is distributed WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef INCTREE_H
#define INCTREE_H

typedef enum {HostNotSet=0, HostSingle, HostAllShared, HostComplex} host_stat_type;

struct Algebra;

#include <stdio.h>
#include "oops/parutil.h"


#define   ON_EVERY_HOST         -1
#define   ID_TREE_NOT_SET       -99


class Tree {
 public: 
  Tree **sons;
  int nb_sons;
  int begin;
  int end;
  int index;
  int nb_nodes;
  bool local;
  Tree *above;
  
  struct Algebra *nodeOfAlg;

#ifdef WITH_MPI
    int id;
    int id_first;
    int id_last;
    MPI_Comm comm;

#ifndef NOOLD
    short *single_first;
    short *single_last;
#endif
    short *share_first;
    short *share_last;
    int   share_sz;
    host_stat_type host_stat;
    MPI_Datatype *bcast_type;

#ifndef NOOLD
    MPI_Datatype reduce_type;
    int bufsize;
#endif
#endif

    Tree(int begin = -1, int end = -1, int nbsons = 0);
    ~Tree();

    void print(FILE *out, const char *name);
    int writeToFile(FILE*);
    static Tree* readFromFile(FILE *f);
    void check();
    bool isIdentical(Tree *T2);

    void setIndex();
    void setLeavesLocal(Tree *Parent);

#ifdef WITH_MPI
    void setProcLocal();
    void setAllNodesBelow(bool local, int id, host_stat_type host_stat);
    void setParBcastTypes();
    void setParBcastTypesNew();
#endif
};

int
comp_nodes(Tree *T1, Tree *T2);

#endif
