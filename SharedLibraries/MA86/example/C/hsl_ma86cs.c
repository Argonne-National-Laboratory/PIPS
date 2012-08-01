/* hsl_ma86cs.c */
#include <stdio.h>
#include <stdlib.h>
#include "hsl_mc69c.h"
#include "hsl_ma86c.h"

/* Simple code to illustrate use of hsl_ma86 */
int main(void) {

   /* Derived types */
   void *keep;
   struct ma86_control control;
   struct ma86_info info;

   int matrix_type, i, n, flag, more;
   int *ptr, *row, *order;
   float complex *val, *x;
   float re, im;

   /* Read the lower triangle of the matrix */
   scanf("%d %d\n", &matrix_type, &n);
   ptr = (int *) malloc(sizeof(int)*(n+1));
   for(i=0; i<n+1; i++) scanf("%d", &ptr[i]);
   row = (int *) malloc(sizeof(int)*(ptr[n]));
   for(i=0; i<ptr[n]; i++) scanf("%d", &row[i]);
   val = (float complex*) malloc(sizeof(float complex)*(ptr[n]));
   for(i=0; i<ptr[n]; i++) {
      scanf(" (%f,%f)", &re, &im);
      val[i] = re + im*I;
   }
   /* Read the right hand side */
   x = (float complex*) malloc(sizeof(float complex)*n);
   for(i=0; i<n; i++) {
      scanf(" (%f,%f)", &re, &im);
      x[i] = re + im*I;
   }

   /* Use the input order */
   order = (int *) malloc(sizeof(int)*n);
   for(i=0; i<n; i++) order[i] = i;

   /* Uncomment the following lines to enable checking (performance overhead) */
   /*
   flag = mc69_verify(6, matrix_type, 0, n, n, ptr, row, NULL, &more);
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
   ma86_factor(matrix_type, n, ptr, row, val, order, &keep, &control, &info, 
      NULL);
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
   for(i=0; i<n; i++) printf("(%5.3f,%5.3f) ", creal(x[i]), cimag(x[i]));
   printf("\n");

   /* Clean up */
   ma86_finalise(&keep, &control);
   free(ptr); free(row); free(val);
   free(order);
   free(x);

   return 0;
}
