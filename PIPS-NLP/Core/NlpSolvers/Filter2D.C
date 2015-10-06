/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/
 
#include "Filter2D.h"
#include "FilterIPMOption.h"
#include <cstddef>
#include <cmath>

Filter2D::Filter2D()
  : filter_ConNorm(0),
    filter_Obj(1e20),
    nextfilter(NULL)
{}

Filter2D::Filter2D(double ConNorm_in, double Obj_in)
{
	filter_ConNorm =  ConNorm_in;
	filter_Obj =  Obj_in;
	nextfilter = NULL;	
}

Filter2D::~Filter2D()
{
  if(nextfilter)
	delete nextfilter;
  nextfilter = NULL;
}


void
Filter2D::Initialize(FilterIPMOption* FilterIPMOpt){
  if(nextfilter)
  {
	delete nextfilter;
    nextfilter = NULL;
  }
  filter_ConNorm = FilterIPMOpt->ConNorm_max;
  filter_Obj = -1e20;
}


bool
Filter2D::WithinFilter(const double curr_ConNorm, double const curr_Obj)
{
  int infilter = false;
  Filter2D *ip;

  for(ip=this;ip;ip=ip->nextfilter)
  {
  	if( curr_ConNorm > ip->filter_ConNorm && curr_Obj > ip->filter_Obj ){
	  infilter = true;
	  break;
	}
  }
  return infilter;
}


void
Filter2D::UpdateFilter(const double curr_ConNorm, const double curr_Obj)
{
  Filter2D *ip = this;

  while(ip->nextfilter!=NULL){
  	ip=ip->nextfilter;
  }
  
  ip->nextfilter = new Filter2D(curr_ConNorm,curr_Obj);
}


