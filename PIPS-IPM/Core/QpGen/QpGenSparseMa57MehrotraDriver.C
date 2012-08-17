/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGenDriver.h"
#include "QpGenSparseMa57.h"
#include "MehrotraSolver.h"
#include <iostream>
#include "omp.h"

int main( int argc, char *argv[] )
{

	MehrotraSolver  * solver = 0;
  QpGenSparseMa57 * qpgen  = 0;

  // Passing nil arguments to keep old compilers that can't handle
  // explicit template instatiation happy.

	
		double startTime = omp_get_wtime();
 
    int result = qpgen_solve( argc, argv, solver, qpgen );
 
    double endTime = omp_get_wtime();
 
    cout << "Duration in seconds="<< (endTime - startTime) << endl;
 
    
	//clock_t startTime = clock();  //Murat added

  
	
	//clock_t endTime = clock();
	//clock_t clockTicksTaken = endTime - startTime;

	//double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;

	//cout<<"SOLVE TIME FOR THE PROBLEM = " << timeInSeconds<< endl ;

  return result;

}
