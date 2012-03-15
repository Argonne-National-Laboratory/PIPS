/* (c) 2009 Marco Colombo, University of Edinburgh.
 *
 * This file is part of SML.
 *
 * SML is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, using version 3 of the License.
 *
 * SML is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 */

/* This is a Cplex driver for the Structured Modelling Language (SML) */

#include <map>
#include <cstdio>
#include <iostream>
#include <limits> // for numeric_limits
#include <ilcplex/cplex.h>
#include "sml-cplex.h"

using namespace std;

class sparsecol {
  public:
   vector<int> row;
   vector<double> val;

   void add_entry(int r, double v) {
      row.push_back(r);
      val.push_back(v);
   }
   int size() const {
      return row.size();
   }
};

int createLP(CPXENVptr env, CPXLPptr lp, ExpandedModelInterface *root);
int getSolution(CPXENVptr env, CPXLPptr lp, ExpandedModelInterface *root);

/** Solve the problem */
int SML_CPLEX_driver(ExpandedModelInterface *root) {

  int status = 0;

  CPXENVptr env;
  CPXLPptr lp;

  // initialize the CPLEX environment
  printf("Calling Cplex...\n");
  env = CPXopenCPLEX(&status);
  if (!env) {

    char errmsg[1024];

    CPXgeterrorstring(env, status, errmsg);
    fprintf(stderr, "Could not open the CPLEX environment.\n%s", errmsg);

    goto TERMINATE;
  }

  // turn on output to the screen
  status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
  if (status) {
    fprintf(stderr, "Failed to turn on screen indicator, error %d.\n", status);
    goto TERMINATE;
  }

  // create the problem
  lp = CPXcreateprob(env, &status, "lp");
  if (!lp) {
    fprintf(stderr, "Failed to create the problem.\n");
    goto TERMINATE;
  }

  // fill the data of the LP problem
  status = createLP(env, lp, root);
  if (status) {
    fprintf(stderr, "Failed to fill the data of the LP problem.\n");
    goto TERMINATE;
  }

  // write the deterministic equivalent in mps format
  if (true) {
    status = CPXwriteprob(env, lp, "smps.lp", NULL);
    if (status) {
      fprintf(stderr, "Failed to write the mps file.\n");
      goto TERMINATE;
    }
  }

  // solve the problem
  status = CPXlpopt(env, lp);
  if (status) {
    fprintf(stderr, "Failed to optimize the LP problem.\n");
    goto TERMINATE;
  }

  // retrieve the solution
  status = getSolution(env, lp, root);
  if (status) {
    fprintf(stderr, "Failed to retrieve the solution.\n");
    goto TERMINATE;
  }

 TERMINATE:

  // free up the problem as allocated by CPXcreateprob, if necessary
  if (lp) {
    status = CPXfreeprob(env, &lp);
    if (status)
      fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
  }

  // free up the CPLEX environment, if necessary
  if (env) {
    status = CPXcloseCPLEX(&env);

    if (status) {

      char errmsg[1024];

      CPXgeterrorstring(env, status, errmsg);
      fprintf(stderr, "Could not close the CPLEX environment.\n%s", errmsg);
    }
  }

  return status;
}

/*
 * Each model block looks like the following:
 * A_1         C_1
 *     A_2     C_2
 *         A_3 C_3
 * B_1 B_2 B_3  D
 *
 * model->getJacobianOfIntersection(model) returns  D
 * model->getJacobianOfIntersection(A_i)   returns B_i
 * A_i->getJacobianOfIntersection(model)   returns C_i
 */
int createLP(CPXENVptr env, CPXLPptr lp, ExpandedModelInterface *root) {

  const double inf = numeric_limits<double>::infinity();

  int status = 0;

  int nRows = 0, nCols = 0, nNonz = 0;

  // count the total number of rows and columns
  for (ExpandedModelInterface::child_iterator i = root->cbegin();
       i != root->cend(); ++i) {

    ExpandedModelInterface *model = *i;
    nRows += model->getNLocalCons();
    nCols += model->getNLocalVars();
  }

  // allocate arrays
  double obj[nCols], rhs[nRows], blo[nCols], bup[nCols];
  char rowType[nRows], **rwnames = NULL, **clnames = NULL;

  map<string,int> row_offset;
  map<string,int> col_offset;
  map<int,sparsecol> cols;
  int next_row = 0, next_col = 0;

  // Setup row and column numbering, together with bounds and constraint types
  for (ExpandedModelInterface::child_iterator i = root->cbegin();
       i != root->cend(); ++i) {

    ExpandedModelInterface *model = *i;
    row_offset[model->getName()] = next_row;
    next_row += model->getNLocalCons();
    col_offset[model->getName()] = next_col;
    next_col += model->getNLocalVars();

    // Name rows and identify constraint types
    {
      int nLocalCons = model->getNLocalCons();
      double lwr_bnds[nLocalCons];
      double upr_bnds[nLocalCons];
      model->getRowBounds(lwr_bnds, upr_bnds);

      int offset = row_offset[model->getName()];
      double *rr = rhs + offset;
      double *lb = lwr_bnds;
      double *ub = upr_bnds;

      for (int j = offset; j < next_row; ++j,++lb,++ub, ++rr) {

        *rr = 0.0;

        if(*lb==-inf && *ub==inf) {
          rowType[j] = 'N';
        }
        else if(*lb==-inf) {
          rowType[j] = 'L';
          if(*ub!=0) *rr = *ub;
        }
        else if(*ub== inf) {
          rowType[j] = 'G';
          if(*lb!=0) *rr = *lb;
        }
        else if(*lb==*ub) {
          rowType[j] = 'E';
          if(*lb!=0) *rr = *lb;
        }
        else {
          rowType[j] = 'L';
          if(*ub!=0) *rr = *ub;
        }
      }
    }

    // set up the objective and bounds
    {
      int nLocalVars = model->getNLocalVars();
      double objv[nLocalVars];
      double lwrb[nLocalVars];
      double uprb[nLocalVars];

      model->getObjGradient(objv);
      model->getColLowBounds(lwrb);
      model->getColUpBounds(uprb);

      int offset = col_offset[model->getName()];
      double *op = obj + offset;
      double *lb = blo + offset;
      double *ub = bup + offset;

      for (int k = 0; k < nLocalVars; ++k, ++op, ++lb, ++ub) {
        *op = objv[k];
        *lb = lwrb[k];
        *ub = uprb[k];
      }
    }
  }

#ifdef DEBUG
  cout << "\nRHS:\n";
  for (int i = 0; i < nRows; ++i) {
    cout << rowType[i] << "  " << rhs[i] << endl;
  }

  cout << "\nOBJ:\n";
  for (int i = 0; i < nCols; ++i) {
    cout << i << "  " << obj[i] << " lo: " << blo[i] << "  up: " << bup[i] <<endl;
  }
#endif

  // create the new columns
  status = CPXnewrows(env, lp, nRows, rhs, rowType, NULL, rwnames);
  if (status)
    goto TERMINATE;

  // add the columns
  status = CPXnewcols(env, lp, nCols, obj, blo, bup, NULL, clnames);
  if (status)
    goto TERMINATE;

  // Now aquire matrix information by ascending and descending the tree
  // Recall that model->getJacobianOfIntersection(jm) returns the
  // block in the current model's row using variables in block jm
  for (ExpandedModelInterface::child_iterator i = root->cbegin();
       i != root->cend(); ++i) {

    ExpandedModelInterface *model = *i;

    // traverse the children
    for (ExpandedModelInterface::child_iterator j = model->cbegin();
         j != model->cend(); ++j) {

      ExpandedModelInterface *jm = *j;

      // skip if ...
      if (model->getNzJacobianOfIntersection(jm) == 0)
        continue;

      int nz = model->getNzJacobianOfIntersection(jm);
      int colbeg[jm->getNLocalVars()+1];
      int collen[jm->getNLocalVars()];
      int rownbs[nz];
      double elts[nz];
      nNonz += nz;

      model->getJacobianOfIntersection(jm, colbeg, collen, rownbs, elts);
      int offset = row_offset[model->getName()];

      for (unsigned int p = 0; p < jm->getLocalVarNames().size(); ++p) {
        sparsecol &col = cols[col_offset[jm->getName()]+p];
        for (int k = colbeg[p]; k < colbeg[p] + collen[p]; ++k) {
          col.add_entry(rownbs[k] + offset, elts[k]);
        }
      }
    }

    // traverse the ancestors
    for (ExpandedModelInterface::ancestor_iterator j = model->abegin();
         j != model->aend(); ++j) {

      ExpandedModelInterface *jm = *j;

      // skip if ...
      if(model->getNzJacobianOfIntersection(jm) == 0)
        continue;

      int nz = model->getNzJacobianOfIntersection(jm);
      int colbeg[jm->getNLocalVars()+1];
      int collen[jm->getNLocalVars()];
      int rownbs[nz];
      double elts[nz];
      nNonz += nz;

      model->getJacobianOfIntersection(jm, colbeg, collen, rownbs, elts);
      int offset = row_offset[model->getName()];

      for (unsigned int p = 0; p < jm->getLocalVarNames().size(); ++p) {
        sparsecol &col = cols[col_offset[jm->getName()]+p];
        for (int k = colbeg[p]; k < colbeg[p] + collen[p]; ++k) {
          col.add_entry(rownbs[k] + offset, elts[k]);
        }
      }
    }
  }

  // convert col into sparse representation for cplex
  {
    int rownmbs[nNonz];
    int colnmbs[nNonz];
    double acoeffs[nNonz];
    int nz = 0;

    for (map<int,sparsecol>::iterator i = cols.begin(); i != cols.end(); ++i) {
      int colidx = i->first;
      sparsecol &col = i->second;
      for (int j = 0; j < col.size(); ++j) {
        rownmbs[nz] = col.row[j];
        colnmbs[nz] = colidx;
        acoeffs[nz] = col.val[j];
        ++nz;
      }
    }

    // create the list of coefficients
    status = CPXchgcoeflist(env, lp, nNonz, rownmbs, colnmbs, acoeffs);
    if (status)
      goto TERMINATE;
  }

 TERMINATE:

  return status;
}

void fillSmlSolution(ExpandedModelInterface *node,
                     double *primal, double *dual,
                     double *slack,  double *rcosts,
                     map<string,int> &row_offset,
                     map<string,int> &col_offset) {

  int row_off = row_offset[node->getName()];
  int col_off = col_offset[node->getName()];

  // pass the solution back to sml
  node->setPrimalSolColumns(primal + col_off);
  node->setDualSolColumns(rcosts + col_off);
  node->setPrimalSolRows(slack + row_off);
  node->setDualSolRows(dual + row_off);

  // recurse down the rest of the tree
  for (int i = 0; i < node->children.size(); ++i) {
    ExpandedModelInterface *model = node->children[i];

    fillSmlSolution(model, primal, dual, slack, rcosts,
                    row_offset, col_offset);
  }
}

/** Retrieve the solution information */
int getSolution(CPXENVptr env, CPXLPptr lp, ExpandedModelInterface *root) {

  int status, solstat, iters;
  double obj;

  int numRows = CPXgetnumrows(env, lp);
  int numCols = CPXgetnumcols(env, lp);

  double *primal = new double[numCols];
  double *dual   = new double[numRows];
  double *slack  = new double[numRows];
  double *rcosts = new double[numCols];

  map<string,int> row_offset;
  map<string,int> col_offset;
  int next_row = 0, next_col = 0;

  if (primal == NULL || dual   == NULL ||
      slack  == NULL || rcosts == NULL) {
    status = CPXERR_NO_MEMORY;
    fprintf(stderr, "Could not allocate memory for the solution.\n");
    goto TERMINATE;
  }

  // retrieve the solution
  status = CPXsolution(env, lp, &solstat, &obj, primal, dual, slack, rcosts);
  if (status)
    goto TERMINATE;

  // retrieve the number of iterations
  iters = CPXgetitcnt(env, lp);

  // output the solution to the screen
  printf("\nObjective function value: %.8f (after %d iterations).\n",
	 obj, iters);

  // Setup row and column numbering, together with bounds and constraint types
  for (ExpandedModelInterface::child_iterator i = root->cbegin();
       i != root->cend(); ++i) {

    ExpandedModelInterface *model = *i;
    row_offset[model->getName()] = next_row;
    next_row += model->getNLocalCons();
    col_offset[model->getName()] = next_col;
    next_col += model->getNLocalVars();
  }

  // pass the solution back to sml
  fillSmlSolution(root, primal, dual, slack, rcosts, row_offset, col_offset);

 TERMINATE:

  // clean up
  delete[] primal;
  delete[] dual;
  delete[] slack;
  delete[] rcosts;

  return status;
}
