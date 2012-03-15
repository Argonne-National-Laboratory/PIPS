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

#ifndef CALLBACK_H
#define CALLBACK_H

typedef struct {

  void *id;
  int max_nz;

  int nz;
  int *row_nbs;
  int *col_beg;
  int *col_len;
  double *element;

} CallBackInterfaceType;


typedef void (*CallBackFunction) (CallBackInterfaceType *cbi);

void 
CallBackVoid(CallBackInterfaceType *cbi);

#endif /* CALLBACK_H */
