/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef VECTORUTILITIES
#define VECTORUTILITIES

#include <iostream>
#include <fstream>
using namespace std;

void set_to_zero( double v[], int n, int stride );
void writef_to_stream( double v[], int n, int stride,
		       ostream& out, const char format[] );
void set_to_constant( double v[], int n, int stride, double c );
void add_constant( double v[], int n, int stride, double c );
double stepbound( double v[], int n, int incv,
				  double s[], int incs, double max );
double stepbound( double v[], int n, int incv,
				  double s[], int incs,
				  double b[], int incb, double u[], int incu,
				  double max );
void axdzpy( int n, double alpha,
	     double x[], int incx, double z[], int incz,
	     double y[], int incy );

double find_blocking( double w[],     int n, int incw,
		    double wstep[],        int incwstep,
		    double u[],            int incu,
		    double ustep[],        int incustep,
		    double maxStep,
		    double *w_elt,          double *wstep_elt,
		    double *u_elt,          double *ustep_elt,
		    int& first_or_second );

void find_blockingPD( double w[],     int n, int incw,
		    double wstep[],        int incwstep,
		    double u[],            int incu,
		    double ustep[],        int incustep,
		    double *w_elt,         double *wstep_elt,
		    double *u_elt,         double *ustep_elt,
		    int& first_or_second,
		    double * alphaPri, double * alphaDual);

void _SetComponentFromMaxXorY( double* z_in, double* x_in, double *y_in, 
				int Start, int End,int xStart, int xEnd, int yStart, int yEnd, double *select = NULL);

void _SetComponentFromMaxXorConstant( double* z_in, double* x_in, double y_in, 
				int Start, int End,int xStart, int xEnd);

void _copyV1intoV2_FromTo(double* V2_in, double* V1_in,int V2Start, int V2End,int V1Start, int V1End);


#endif
