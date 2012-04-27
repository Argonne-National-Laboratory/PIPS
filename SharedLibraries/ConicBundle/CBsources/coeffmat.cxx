/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/coeffmat.cxx

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

/* ****************************************************************************

    Copyright (C) 2004-2010  Christoph Helmberg

    ConicBundle, Version 0.3.8
    File:  CBsources/coeffmat.cxx

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



#include <string.h>
#include <stdlib.h>
#include <cctype>
#include "coeffmat.hxx"
#include "cmsymden.hxx"
#include "cmsymspa.hxx"
#include "cmgramde.hxx"
#include "cmgramsp.hxx"
#include "cmrankdd.hxx"
#include "cmranksd.hxx"
#include "cmrankss.hxx"
#include "cmsingle.hxx"

 

namespace ConicBundle {

Coeffmat* coeffmat_read(std::istream& in)
{
 char name[80];
 in>>std::ws;
 if (!in) {
     if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<"*** ERROR: coeffmat_read(): input stream broken"<<std::endl;
     return 0;
 }
 int cnt=0;
 while((in)&&(!isspace(in.peek()))&&(cnt<80)){
   in.get(name[cnt++]);
 }
 if ((!in)||(cnt==80)){
     if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<"*** ERROR: coeffmat_read(): failed in reading name of constraint"<<std::endl;
     return 0;
 }
 name[cnt]=0;
 Coeffmat* p=0;
 if (!strcmp(name,"SYMMETRIC_DENSE")) p=new CMsymdense(in);
 else if (!strcmp(name,"SYMMETRIC_SPARSE")) p=new CMsymsparse(in);
 else if (!strcmp(name,"GRAM_DENSE")) p=new CMgramdense(in);
 else if (!strcmp(name,"GRAM_SPARSE")) p=new CMgramsparse(in);
 else if (!strcmp(name,"LOWRANK_DENSE_DENSE")) p=new CMlowrankdd(in);
 else if (!strcmp(name,"LOWRANK_SPARSE_DENSE")) p=new CMlowranksd(in);
 else if (!strcmp(name,"LOWRANK_SPARSE_SPARSE")) p=new CMlowrankss(in);
 else if (!strcmp(name,"SINGLETON")) p=new CMsingleton(in);
 else {
     if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<"*** ERROR: coeffmat_read(): unknown constraint name :"<<name<<std::endl;
     in.clear(std::ios::failbit);
     return 0;
 }
 if (p==0) {
     if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<"*** ERROR: coeffmat_read():";
     if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<" failed in reading a constraint of type "<<name<<std::endl;
     in.clear(std::ios::failbit);
     return 0;
 }
 return p;
}

}

