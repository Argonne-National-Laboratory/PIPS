/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/CB_CSolver.hxx

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



#ifndef CONICBUNDLE_CB_CSOLVER_HXX
#define CONICBUNDLE_CB_CSOLVER_HXX

#include <map>
#include "cb_cinterface.h"
#include "MatBSolver.hxx"
#include "CFunction.hxx"


class CB_CSolver
{
public:
  std::map<void*, ConicBundle::CFunction*> funmap;
  ConicBundle::MatrixBSolver* solver;
 
  CB_CSolver(bool no_bundle);

  ~CB_CSolver();
	  
};

#endif

