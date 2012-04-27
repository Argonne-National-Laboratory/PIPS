/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  c_main.c

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



#include <stdio.h>
#include "cb_cinterface.h"

/* prototypes for two oracle functions */

int eval_fun0( /* input: */
	      void* function_key,
	      double *dual, /* argument/Lagrange multipliers */
	      double relprec,
	      int max_new_subg,
	      
	      /* output: */ 
	      double *objective_value,
	      int* new_subg,
	      double *subgval,
	      double *subgradient,
	      double *primal           /* may be zero pointer */
	      );


int eval_fun1( /* input: */
	      void* function_key,
	      double *dual, /* argument/Lagrange multipliers */
	      double relprec,
	      int max_new_subg,
	      
	      /* output: */ 
	      double *objective_value,
	      int* new_subg,
	      double *subgval,
	      double *subgradient,
	      double *primal           /* may be zero pointer */
	      );


/* implementation of the two functions */

int eval_fun0( /* input: */
	      void* function_key,
	      double *dual, /* argument/Lagrange multipliers */
	      double relprec,
	      int max_new_subg,
	      
	      /* output: */ 
	      double *objective_value,
	      int* new_subg,
	      double *subgval,
	      double *subgradient,
	      double *primal           /* may be zero pointer */
	      )
{
  double s0,s1;
  s0=dual[0]-100.;
  subgradient[0]=2.*s0;
  s1=dual[1]-100.;
  subgradient[1]=2.*s1;
  if (primal!=0){
    primal[0]=subgradient[0];
    primal[1]=subgradient[1];
  }
  *objective_value=s0*s0+s1*s1;
  subgval[0]=*objective_value;
  *new_subg=1;
  if (max_new_subg>1){
    /* two times the same subgradient for checking the bundle update */
    subgradient[2]=subgradient[0];
    subgradient[3]=subgradient[1];
    subgval[1]=subgval[0];
    *new_subg=2;
    if (primal!=0){
      primal[2]=subgradient[0];
      primal[3]=subgradient[1];
    }
  }
    
  return 0;
}
  
int eval_fun1( /* input: */
	      void* function_key ,
	      double *dual, /* argument/Lagrange multipliers */
	      double relprec ,
	      int max_new_subg,
	      
	      /* output: */ 
	      double *objective_value,
	      int* new_subg,
	      double *subgval,
	      double *subgradient,
	      double *primal           /* may be zero pointer */
	      )
/* (y0-100)^2+(y1-98)^2 */
{
  double s0,s1;
  int i;
  s0=dual[0]-100.;
  subgradient[0]=2.*s0;
  s1=dual[1]-98;
  subgradient[1]=2.*s1;
  *objective_value=s0*s0+s1*s1;
  *subgval=*objective_value;
  *new_subg=1;
  if (max_new_subg>1){
    double d0=0.01;
    double d1=0.01;
    s0=dual[0]+d0-100.;
    subgradient[2]=2.*s0;
    s1=dual[1]+d1-98;
    subgradient[3]=2.*s1;
    subgval[1]=s0*s0+s1*s1-subgradient[2]*d0-subgradient[3]*d1;
    *new_subg=2;
  }
  if (max_new_subg>2){
    double d0=0.01;
    double d1=-0.01;
    s0=dual[0]+d0-100.;
    subgradient[4]=2.*s0;
    s1=dual[1]+d1-98;
    subgradient[5]=2.*s1;
    subgval[2]=s0*s0+s1*s1-subgradient[4]*d0-subgradient[5]*d1;
    *new_subg=3;
  }
  if (max_new_subg>3){
    double d0=-0.01;
    double d1=0.01;
    s0=dual[0]+d0-100.;
    subgradient[6]=2.*s0;
    s1=dual[1]+d1-98;
    subgradient[7]=2.*s1;
    subgval[3]=s0*s0+s1*s1-subgradient[6]*d0-subgradient[7]*d1;
    *new_subg=4;
  }
  if (max_new_subg>4){
    double d0=-0.01;
    double d1=-0.01;
    s0=dual[0]+d0-100.;
    subgradient[8]=2.*s0;
    s1=dual[1]+d1-98;
    subgradient[9]=2.*s1;
    subgval[4]=s0*s0+s1*s1-subgradient[8]*d0-subgradient[9]*d1;
    *new_subg=5;
  }
  if (primal!=0) {
    for (i=0;i<*new_subg;i++){
      primal[2*i]=subgradient[2*i];
      primal[2*i+1]=subgradient[2*i+1];
    }
  }
  return 0;
}
  



int main(void)
{
  cb_problemp p;
  double y[2];
  double x0[2];
  double x1[2];
  double obj;
  double u;
  int cnt;


  p=cb_construct_problem(0);
  if (p==0){
    printf("construct_problem failed\n");
    return 1;
  }
  if (cb_init_problem(p,2,0,0)){
    printf("init_problem failed\n");
    return 1;
  }
  if (cb_add_function(p,(void *)eval_fun0,eval_fun0,0,2)){
    printf("add eval_fun0 failed\n");
    return 1;
  }
  if (cb_add_function(p,(void *)eval_fun1,eval_fun1,0,2)){
    printf("add eval_fun1 failed\n");
    return 1;
  }
    
  cb_set_print_level(p,1);
  cb_set_term_relprec(p,1e-5);
  cb_set_max_bundlesize(p,(void *)eval_fun0,10);
  cb_set_max_bundlesize(p,(void *)eval_fun1,10);
  cb_set_max_new_subgradients(p,(void *)eval_fun0,2);
  cb_set_max_new_subgradients(p,(void *)eval_fun1,5);
  /* cb_set_next_weight(p,5.); */
  
  cnt=1;

  do {

    int retval;

    
    /* make a descent step */
    if ((retval=cb_do_descent_step(p))){
      printf("cb_do_descent_step returned %d\n",retval);
      return 1;
    }

    
    /* get solution information */
    obj=cb_get_objval(p);
    if ((retval=cb_get_center(p,y))){
      printf("cb_get_center returned %d\n",retval);
      return 1;
    }
    u=cb_get_last_weight(p);

    printf("%d: %f [%f,%f] %f",cnt,obj,y[0],y[1],u);


    /* get primal information */
    if ((retval=cb_get_approximate_primal(p,(void*)eval_fun0,x0))){
      printf("cb_get_primal returned %d for function 0\n",retval);
      return 1;
    }
    printf("\n xf0: [%f,%f]",x0[0],x0[1]);

    if ((retval=cb_get_approximate_primal(p,(void*)eval_fun1,x1))){
      printf("cb_get_primal returned %d for function 1\n",retval);
      return 1;
    }      
    printf("\n xf1: [%f,%f]",x1[0],x1[1]);

    /* set some paramters */
    /* cb_set_min_weight(p,10.); */
    /* cb_set_max_weight(p,10.); */

    printf("\n");
    cnt++;

  } while (!cb_termination_code(p));

  cb_print_termination_code(p);

  cb_destruct_problem(&p);

  return 0;
}
