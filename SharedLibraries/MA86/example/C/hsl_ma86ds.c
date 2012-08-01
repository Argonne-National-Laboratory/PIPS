/* hsl_ma86ds.c */
#include <stdio.h>
#include <stdlib.h>
#include "hsl_mc69d.h"
#include "hsl_ma86d.h"

/* Simple code to illustrate use of hsl_ma86 */
int main(void) {
   void *keep;
   struct ma86_control control;
   struct ma86_info info;

   int i, n, flag, more;
   int *ptr, *row, *order;
   double *val, *x;

   /* Read the lower triangle of the matrix */
   scanf("%d\n", &n);
   ptr = (int *) malloc(sizeof(int)*(n+1));
   for(i=0; i<n+1; i++) scanf("%d", &ptr[i]);
   row = (int *) malloc(sizeof(int)*(ptr[n]));
   for(i=0; i<ptr[n]; i++) scanf("%d", &row[i]);
   val = (double *) malloc(sizeof(double)*(ptr[n]));
   for(i=0; i<ptr[n]; i++) scanf("%lf", &val[i]);
   /* Read the right hand side */
   x = (double *) malloc(sizeof(double)*n);
   for(i=0; i<n; i++) scanf("%lf", &x[i]);

   /* Use the input order */
   order = (int *) malloc(sizeof(int)*n);
   for(i=0; i<n; i++) order[i] = i;

   /* Uncomment the following lines to enable checking (performance overhead) */
   /*
   flag = mc69_verify(6, HSL_MATRIX_REAL_SYM_INDEF, 0, n, n, ptr, row, NULL,
      &more);
   if(flag != 0) {
      printf("Matrix not in HSL standard format. flag, more = %i, %i\n", 
         flag, more);
      return 1;
   }
   */

   /* Set up control */
   ma86_default_control(&control);

   /* Analyse */
   ma86_analyse(n, ptr, row, order, &keep, &control, &info);
   if(info.flag < 0) {
      printf("Failure during analyse with info.flag = %i\n", info.flag);
      return 1;
   }

   /* Factor */
   ma86_factor(n, ptr, row, val, order, &keep, &control, &info, NULL);
   if(info.flag < 0) {
      printf("Failure during factor with info.flag = %i\n", info.flag);
      return 1;
   }

   /* Solve */
   ma86_solve(0, 1, n, x, order, &keep, &control, &info, NULL);
   if(info.flag < 0) {
      printf("Failure during solve with info.flag = %i\n", info.flag);
      return 1;
   }

   printf(" Computed solution:\n");
   for(i=0; i<n; i++) printf("%10.3lf ", x[i]);
   printf("\n");

   /* Clean up */
   ma86_finalise(&keep, &control);
   free(ptr); free(row); free(val);
   free(order);
   free(x);

   return 0;
}
