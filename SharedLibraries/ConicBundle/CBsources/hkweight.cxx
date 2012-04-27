/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/hkweight.cxx

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



#include "hkweight.hxx"
#include "mymath.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

BundleHKweight::BundleHKweight(Real mRin)
{
  out=0;
  print_level=1;
  mR=mRin;
  clear();
}

void BundleHKweight::set_defaults()
{
  mR=0.5;
} 

void BundleHKweight::clear()
{
  modelmax=CB_minus_infinity;
  minweight=maxweight=weight=-1.;
  weightchanged=0; 
  iweight=0;
  epsweight=1e30;
  next_weight_set=0;
}
  
int BundleHKweight::init(const Matrix& subgrad,Real normsubg2)
{
 modelmax=CB_minus_infinity;
 if (!next_weight_set) {
   if (subgrad.rowdim()==0)
     weight=1.;
   else {
     Real d=sqrt(normsubg2);
     if (d<subgrad.rowdim()*1e-10)
       weight=1.;
     else
       weight=max(d,1e-4);
   }
 }
 if (minweight<=0) minweight=max(1e-10*weight,1e-10);
 weight=max(minweight,weight);
 if (maxweight>0) weight=min(weight,maxweight);  
 epsweight=1e30;
 iweight=0;
 weightchanged=1;
 return 0;
}

int BundleHKweight::init(const Matrix& subgrad)
{
 modelmax=CB_minus_infinity;
 if (!next_weight_set){
   if (subgrad.rowdim()==0)
     weight=1.;
   else {
     Real d=norm2(subgrad);
     if (d<subgrad.rowdim()*1e-10)
       weight=1.;
     else
       weight=max(norm2(subgrad),1e-4);
   }
 }
 if (minweight<=0) minweight=max(1e-10*weight,1e-10);
 weight=max(minweight,weight);
 if (maxweight>0) weight=min(weight,maxweight);  
 epsweight=1e30;
 iweight=0;
 weightchanged=1;
 return 0;
}

int BundleHKweight::prob_changed(const Matrix& /* subgrad */,Real normsubg2)
{
 modelmax=CB_minus_infinity;
 if (!next_weight_set) weight=max(weight,sqrt(normsubg2));
 if (minweight<=0) minweight=max(1e-10*weight,1e-10);
 weight=max(minweight,weight);
 if (maxweight>0) weight=min(weight,maxweight);  
 weightchanged=1;
 /*
 weight=sqrt(normsubg2);
 if (minweight<=0) minweight=1e-10*weight;
 else weight=max(minweight,weight);
 if (maxweight>0) weight=min(weight,maxweight);  
 weightchanged=1;
 */
 epsweight=1e30;
 iweight=0;
 return 0;
}

Real BundleHKweight::get_weight() const
{ return weight; }
    
int BundleHKweight::weight_changed() const
{ return weightchanged; }


int BundleHKweight::descent_update(Real newval,Real oldval,Real modelval,
				   const Matrix& /* y */, const Matrix& /* newy */, Real /* norm2subg */)
{
 if (weight<0) return 1;
 next_weight_set=0;

 modelmax=max(modelmax,modelval);
 Real oldweight=weight;
 Real weightint=2*weight*(1.-(oldval-newval)/(oldval-modelval));
 if ((out)&&(print_level>0)) (*out)<<"  serious step, i_u="<<iweight<<std::flush;
 if (((oldval-newval)>mR*(oldval-modelval))&&(iweight>0)){
     if ((out)&&(print_level>0)) (*out)<<" uint="<<weightint<<std::flush;
     weight=weightint;
 }
 else if (iweight>3){ //there were four, now is the 5th consecutive serious 
     if ((out)&&(print_level>0)) (*out)<<" i_u>3 u/2 "<<std::flush;
     weight/=2.;
 }
 else if (newval<modelmax){
     if ((out)&&(print_level>0)) (*out)<<" nv<linmax u/2 "<<std::flush;
     weight/=2.;
 }
 weight=max(oldweight/10.,weight);
 if (minweight>0) weight=max(weight,minweight);
 if ((out)&&(print_level>0)) (*out)<<" unew="<<weight<<std::endl;
 epsweight=max(epsweight,2*(oldval-modelval));
 iweight=max(iweight+1,Integer(1));
 if (weight<oldweight) {
     weightchanged=1;
     iweight=1;
     modelmax=CB_minus_infinity;
 }
 else {
     weightchanged=0;
 }
 return 0;
}


int BundleHKweight::nullstep_update(
		       Real newval,Real oldval,Real modelval,Real lin_approx,
                       const Matrix& /* y */, const Matrix& /* newy */, Real normsubg2)
{
 if (weight<0) return 1;
 next_weight_set=0;
 
 Real oldweight=weight;
 Real weightint=2*weight*(1.-(oldval-newval)/(oldval-modelval));

 epsweight=min(epsweight,sqrt(normsubg2*weight)+oldval-modelval-normsubg2); 
 if ((out)&&(print_level>0)) (*out)<<"  null step, i_u="<<iweight<<" eps_v="<<epsweight<<" lapprox="<<lin_approx<<std::flush;
 
 if ((lin_approx>max(epsweight,10.*(oldval-modelval)))&&(iweight<-3)){
     if ((out)&&(print_level>0)) (*out)<<" uint="<<weightint<<std::flush;
     weight=weightint;
 }
 weight=min(weight,10.*oldweight);   
 if (maxweight>0) weight=min(weight,maxweight);
 if ((out)&&(print_level>0)) (*out)<<" unew="<<weight<<std::endl;
 iweight=min(iweight-1,Integer(-1));
 if (weight>oldweight) {
     iweight=-1;
     weightchanged=1;
 }
 else {
     weightchanged=0;
 }
 return 0;
}

}

