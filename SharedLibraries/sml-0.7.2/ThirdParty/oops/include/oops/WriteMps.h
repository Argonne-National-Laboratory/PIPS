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
#ifndef WRITEMPS_H
#define WRITEMPS_H

#include "oops/Vector.h"
#include "oops/Algebra.h"

int
Write_MpsFile(FILE *out, Algebra * AlgAug, Vector *rhs, Vector *obj,
              Vector *upb, Vector *lob,
	      int classic_mps, char **colnames, char **rownames);

int
WriteMps(const char *filename, Algebra *AlgAug,
	 Vector *rhs, Vector *obj, Vector *upb);

#endif /* WRITEMPS_H */
