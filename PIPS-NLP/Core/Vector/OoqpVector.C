/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "OoqpVector.h"
#include <cmath>

OoqpVector::OoqpVector( int n_ )
{
  n = n_;
}

OoqpVector::~OoqpVector()
{
}


// default is inf norm
double 
OoqpVector::Norm(const int Norm_Type, OoqpVector *vec2)
{
  double wrkDouble1,wrkDouble2;
  
  if(vec2 == NULL){
    switch(Norm_Type){
  	  case PIPS_NORM_INF:
		return this->infnorm();
	  case PIPS_NORM_ONE:
		return this->onenorm();
	  case PIPS_NORM_TWO:
		return this->twonorm();
	  default: 
		return this->infnorm();
    }
  }else{
    switch(Norm_Type){
  	  case PIPS_NORM_INF:
	  	wrkDouble1 = this->infnorm();
		wrkDouble2 = vec2->infnorm();
		return (wrkDouble1>wrkDouble2)?wrkDouble1:wrkDouble2;
	  case PIPS_NORM_ONE:
		wrkDouble1 = this->onenorm();
		wrkDouble2 = vec2->onenorm();
		return wrkDouble1+wrkDouble2;
	  case PIPS_NORM_TWO:
		wrkDouble1 = this->twonorm(); wrkDouble1*=wrkDouble1;
		wrkDouble2 = vec2->twonorm(); wrkDouble2*=wrkDouble2;
		return sqrt(wrkDouble1+wrkDouble2);
	  default: 
		wrkDouble1 = this->infnorm();
		wrkDouble2 = vec2->infnorm();
		return (wrkDouble1>wrkDouble2)?wrkDouble1:wrkDouble2;

    }
  }
}


double 
OoqpVector::Norm(const int Norm_Type, const double norm1, const double norm2)
{
  double wrkDouble1,wrkDouble2;
  
  switch(Norm_Type){
    case PIPS_NORM_INF:
	  return (norm1>norm2)?norm1:norm2;
	case PIPS_NORM_ONE:
	  return norm1+norm2;
	case PIPS_NORM_TWO:
	  return sqrt(norm1*norm1+norm2*norm2);
	default: 
	  return (norm1>norm2)?norm1:norm2;
  }

}

void absVal(OoqpVector *vec_in){assert("Not Implemented" && 0);};

