/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/qp_solver.hxx

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


#ifndef CONICBUNDLE_QP_SOLVER_HXX
#define CONICBUNDLE_QP_SOLVER_HXX

#include <vector>
#include "qp_block.hxx"

namespace ConicBundle {

class QP_Solver
{
 private:
  //--- problem desription
  CH_Matrix_Classes::Symmatrix Q;  //quadratic cost matrix (positive definite)
  CH_Matrix_Classes::Matrix c;     //linear cost matrix
  CH_Matrix_Classes::Real offset;  //constant added to objective
  std::vector<QP_Block *> blockp; //holds pointers to the blocks 
                //(these are owned by someone else)
                //in principle these have to be set only once,
                //it is important that they are in the same sequence as
                //used for constructing Q so that the variable indices correspond
                //the blocks are not deleted on exit

  CH_Matrix_Classes::Matrix A;     //a (full) blockmatrix, it is collected by calling blockp[i]
  CH_Matrix_Classes::Matrix b;     //primal right hand side, it is collected by calling blockp[i]
 
  //--- termination parameters
  CH_Matrix_Classes::Real lowerbound;
  CH_Matrix_Classes::Real upperbound;
  CH_Matrix_Classes::Real termeps;
  CH_Matrix_Classes::Integer maxiter;

  //--- global variables
  CH_Matrix_Classes::Matrix x,y,z;         //current point

  CH_Matrix_Classes::Real primalval;       //primal objective value
  CH_Matrix_Classes::Real dualval;         //dual objective value

  CH_Matrix_Classes::Real mu;              //barrier parameter
  CH_Matrix_Classes::Matrix Qx;            //holds Q*x for current x

  CH_Matrix_Classes::Integer iter;         //counts number o iterations
  CH_Matrix_Classes::Integer status;       //termination status of last call to solve or update

  //--- temporary variables, global only for memory managment purposes
  CH_Matrix_Classes::Matrix dx,dy,dz;      //step direction
  CH_Matrix_Classes::Symmatrix LDL_Q;      // L*L^T factorization of Q+blockdiag
  CH_Matrix_Classes::Matrix LinvAt;        //=L^-1*A^T
  CH_Matrix_Classes::Symmatrix sysdy;      //system matrix for dy
  CH_Matrix_Classes::Matrix rd;            //dual slack rd=c-Qx-At*y (=-z if feasible)
  CH_Matrix_Classes::Matrix rhs;           //rhs for the system
  CH_Matrix_Classes::Matrix xcorr;
  CH_Matrix_Classes::Matrix tmpvec;        //temporary vector, could be avoided
   
  //--- output
  std::ostream* out;
  int print_level;
  
  //--- private routines
  int iterate();

 public:
  void clear();
  void set_defaults(); 

  QP_Solver(){out=0;print_level=1;clear();set_defaults();}
  ~QP_Solver(){}

  void init_size(CH_Matrix_Classes::Integer maxdim){Q.newsize(maxdim); x.newsize(maxdim,1);}
  void clear_blocks(){blockp.clear();}
  
  const CH_Matrix_Classes::Symmatrix& get_Q(void) const {return Q;}
  const CH_Matrix_Classes::Matrix& get_c(void) const {return c;}
  CH_Matrix_Classes::Real get_offset(void) const {return offset;}
  void add_block(QP_Block* p){blockp.push_back(p);}
  
  void set_termbounds(CH_Matrix_Classes::Real lb,CH_Matrix_Classes::Real ub){lowerbound=lb;upperbound=ub;}
  void set_termeps(CH_Matrix_Classes::Real te){termeps=te;}
  void set_maxiter(CH_Matrix_Classes::Integer mi){maxiter=mi;}
  void set_out(std::ostream* o=0,int pril=1){out=o;print_level=pril;}
  
  CH_Matrix_Classes::Integer get_iter() const {return iter;}
  CH_Matrix_Classes::Integer get_status() const {return status;}
  CH_Matrix_Classes::Real get_termeps() const {return termeps;}
  CH_Matrix_Classes::Integer get_maxiter() const {return maxiter;}
  
  CH_Matrix_Classes::Real get_primalval()const {return primalval;}
  CH_Matrix_Classes::Real get_dualval()const {return dualval; }
  
  int solve(const CH_Matrix_Classes::Symmatrix& Q,const CH_Matrix_Classes::Matrix& c,CH_Matrix_Classes::Real offset);
  // returns status
  //   0 ... if solve terminated correctly for the standard stopping criterion
  //   1 ... if terminated by maxiter or precision tolerances
  //   2 ... if the factorization of the system failed
  //   3 ... if one of the blocks failed 
  
  int resolve();
  // returns status
  //   0 ... if solve terminated correctly for the standard stopping criterion
  //   1 ... if terminated by maxiter or precision tolerances
  //   2 ... if the factorization of the system failed
  //   3 ... if one of the blocks failed 
  
  int update(const CH_Matrix_Classes::Matrix& dc,CH_Matrix_Classes::Real doffset);
  // solves the system for c+=dc offset+=doffset and returns status
  //   0 ... if solve terminated correctly for the standard stopping criterion
  //   1 ... if terminated by maxiter or precision tolerances
  //   2 ... if the factorization of the system failed
  //   3 ... if one of the blocks failed 
  
  std::ostream& save(std::ostream& out) const;
  std::istream& restore(std::istream& in);
};

}

#endif

