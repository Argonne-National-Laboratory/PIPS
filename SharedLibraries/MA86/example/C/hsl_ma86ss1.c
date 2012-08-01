/* hsl_ma86ss1.c */
#include <stdio.h>
#include <stdlib.h>
#include "hsl_mc68i.h"
#include "hsl_mc69s.h"
#include "hsl_ma86s.h"

void stop_on_bad_flag(const char *context, const int flag);

/* Code to illustrate use of more complex hsl_ma86 features */
int main(void) {
   typedef float pkgtype;

   struct mc68_control control68;
   struct mc68_info info68;

   void *keep;
   struct ma86_control control;
   struct ma86_info info;

   int i, j, n, ne, nrhs, lmap, lrow, noor, ndup, flag;
   int *crow, *ccol, *ptr, *row, *order, *map;
   pkgtype *cval, *val, *x, *x2;

   /* Read the first matrix in coordinate format */
   scanf("%d %d\n", &n, &ne);
   ccol = (int *) malloc(sizeof(int)*ne);
   for(i=0; i<ne; i++) scanf("%d", &ccol[i]);
   crow = (int *) malloc(sizeof(int)*ne);
   for(i=0; i<ne; i++) scanf("%d", &crow[i]);
   cval = (pkgtype *) malloc(sizeof(pkgtype)*ne);
   for(i=0; i<ne; i++) scanf("%f", &cval[i]);
   /* Read the first right hand side */
   x = (pkgtype *) malloc(sizeof(pkgtype)*n);
   for(i=0; i<n; i++) scanf("%f", &x[i]);

   /* Convert to HSL standard format */
   ptr = (int *) malloc(sizeof(int)*(n+1));
   lrow = ne; /* maximum size of output matrix cannot exceed size of input */
   lmap = 2*ne; /* large enough for worst case */
   row = (int *) malloc(sizeof(int)*lrow);
   val = (pkgtype *) malloc(sizeof(pkgtype)*lrow);
   map = (int *) malloc(sizeof(int)*lmap);
   flag = mc69_coord_convert(6, HSL_MATRIX_REAL_SYM_INDEF, 0, n, n, ne,
      crow, ccol, cval, ptr, lrow, row, val, &noor, &ndup, &lmap, map);
   stop_on_bad_flag("mc69_coord_convert", flag);

   /* Call mc68 to find a fill reducing ordering (1=AMD) */
   order = (int *) malloc(sizeof(int)*n);
   mc68_default_control(&control68);
   mc68_order(1, n, ptr, row, order, &control68, &info68);
   stop_on_bad_flag("mc68_order", info68.flag);

   /* Initialize control parameters */
   ma86_default_control(&control);

   /* Analyse */
   ma86_analyse(n, ptr, row, order, &keep, &control, &info);
   stop_on_bad_flag("analyse", info.flag);

   /* Factor */
   ma86_factor(n, ptr, row, val, order, &keep, &control, &info, NULL);
   stop_on_bad_flag("factor", info.flag);

   /* Solve */
   ma86_solve(0, 1, n, x, order, &keep, &control, &info, NULL);
   stop_on_bad_flag("solve", info.flag);

   printf(" Computed solution:\n");
   for(i=0; i<n; i++) printf(" %f", x[i]);

   /* Read the values of the second matrix and the new right hand sides */
   for(i=0; i<ne; i++) scanf("%f", &cval[i]);
   scanf("%d\n", &nrhs);
   x2 = (pkgtype *) malloc(sizeof(pkgtype)*n*nrhs);
   for(i=0; i<nrhs; i++) {
      for(j=0; j<n; j++) {
         scanf("%f", &x2[i*n+j]);
      }
   }
   printf("\n");

   /* Convert the values to HSL standard form */
   mc69_set_values(HSL_MATRIX_REAL_SYM_INDEF, lmap, map, cval,
      ptr[n], val);

   /* Perform second factorization and solve */
   ma86_factor_solve(n, ptr, row, val, order, &keep, &control, &info,
      nrhs, n, x2, NULL);
   stop_on_bad_flag("factor_solve", info.flag);

   printf(" Computed solutions:\n");
   for(i=0; i<nrhs; i++) {
      for(j=0; j<n; j++) printf(" %f", x2[i*n+j]);
      printf("\n");
   }

   /* Clean up */
   ma86_finalise(&keep, &control);
   free(crow); free(ccol); free(cval);
   free(ptr); free(row); free(val);
   free(map); free(order);
   free(x); free(x2);

   return 0;
}

void stop_on_bad_flag(const char *context, const int flag) {
   if(0==flag) return;
   printf("Failure during %s  with flag = %d\n", context, flag);
   exit(1);
}
