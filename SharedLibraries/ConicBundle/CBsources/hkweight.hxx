/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/hkweight.hxx

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



#ifndef CONICBUNDLE_HKWEIGHT_HXX
#define CONICBUNDLE_HKWEIGHT_HXX

#include "bundle.hxx"

namespace ConicBundle {

class BundleHKweight: public BundleWeight
{
  CH_Matrix_Classes::Integer iweight;
  CH_Matrix_Classes::Real weight;
  CH_Matrix_Classes::Real epsweight;        //helps in choosing tau
  CH_Matrix_Classes::Real minweight;
  CH_Matrix_Classes::Real maxweight;
  CH_Matrix_Classes::Real modelmax;
  
  int weightchanged;    //1 if last choose_* call modified tau
  int next_weight_set;  //1 if set_nextweight was just called
                        //reset at *_update
  
  CH_Matrix_Classes::Real mR;            //parameter for reduction criterion in serious
  
  std::ostream* out;
  int print_level;
  
public:
  BundleHKweight(CH_Matrix_Classes::Real mRin=.5);
  ~BundleHKweight(){}

  virtual void set_defaults(); 
  //set default values for 'constant' parameters, e.g. minweight and maxweight

  virtual void clear();        
  //reset all adaptive variables and parameters
  
  int init(const CH_Matrix_Classes::Matrix& subgrad,CH_Matrix_Classes::Real normsubg2);
  //compute first tau and set some parameters

  int init(const CH_Matrix_Classes::Matrix& subgrad);
  //compute first tau and set some parameters
  
  int prob_changed(const CH_Matrix_Classes::Matrix& subgrad,CH_Matrix_Classes::Real normsubg2);
  //called if variables were added or deleted
  
  void set_next_weight(CH_Matrix_Classes::Real u)
    //<=0 leaves everything unchanged and does nothing
  { if (u<=0.) return; 
  weight=CH_Matrix_Classes::max(u,1e-10); weightchanged=1;next_weight_set=1;
  modelmax=CB_minus_infinity;iweight=0;}
  
  void set_minweight(CH_Matrix_Classes::Real mw)
    //<=0 means no bound 
  { 
    minweight=mw; 
    if (minweight>0){
      if (weight<minweight) 
	weight=minweight;
      if ((maxweight>0)&&(maxweight<minweight)) 
	maxweight=minweight;
    }
  }
  
  virtual void set_maxweight(CH_Matrix_Classes::Real mw)
    //<=0 means no bound
  { 
    maxweight=mw; 
    if (maxweight>0){
      if (weight>maxweight){
	weight=maxweight;
      }
      if ((minweight>0)&&(minweight>maxweight)) 
	minweight=maxweight;
    }
  }
  
    CH_Matrix_Classes::Real get_weight() const;
  //returns current value of tau
  
  int weight_changed() const;
  //returns 1 if last call of *_update modified current value of tau, else 0
  
  int descent_update(CH_Matrix_Classes::Real newval,CH_Matrix_Classes::Real oldval,CH_Matrix_Classes::Real modelval,
		     const CH_Matrix_Classes::Matrix& y, const CH_Matrix_Classes::Matrix& newy,CH_Matrix_Classes::Real normsubg2);
  //determine next weight after a descent step
  
  int nullstep_update(CH_Matrix_Classes::Real newval,CH_Matrix_Classes::Real oldval,CH_Matrix_Classes::Real modelval,CH_Matrix_Classes::Real lin_approx,
		      const CH_Matrix_Classes::Matrix& y, const CH_Matrix_Classes::Matrix& newy,CH_Matrix_Classes::Real normsubg2);
  //determine next weight after a null step
    
  void set_out(std::ostream* inout=0,int pril=1){out=inout;print_level=pril;}
  
};

}

#endif

