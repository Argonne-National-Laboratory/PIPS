/* (c) 2008,2009 Jonathan Hogg and Andreas Grothey, University of Edinburgh
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
/* This is an MPS driver for the Structured Modelling Language (SML) */

#include <map>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits> // for numeric_limits
#include "sml-mps.h"

using namespace std;

class sparse_col {
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

class mps_bound {
  public:
   string type;
   int id;
   double value;

   mps_bound(const char type_[3], int id_, double value_) :
      type(type_), id(id_), value(value_) {}
};

void writeMps(ostream &out, string name, map<int,char> rows,
   map<int,double> rhs, map<int,double> ranges, list<mps_bound> bnds,
   map<int,sparse_col> cols);

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
void SML_MPS_driver(ExpandedModelInterface *root, string filename) {
   const double inf = numeric_limits<double>::infinity();
   const string obj_name = "obj";

   map<int,char> rows;
   map<int,double> rhs;
   map<int,sparse_col> cols;
   map<int,double> ranges;
   list<mps_bound> bnds;
   int next_row = 1;
   int next_col = 0;

   // Initialise objective
   rows[0] = 'N';

   // Setup row and column numbering, together with bounds and constraint types
   map<string,int> row_offset;
   map<string,int> col_offset;
   for(ExpandedModelInterface::child_iterator i=root->cbegin(); i!=root->cend(); ++i) {
      ExpandedModelInterface *model = *i;
      row_offset[model->getName()] = next_row;
      next_row += model->getNLocalCons();
      col_offset[model->getName()] = next_col;
      next_col += model->getNLocalVars();


      // Name columns, add objective entries and bounds
      {
         double obj[model->getNLocalVars()];
         double lwr_bnds[model->getNLocalVars()];
         double upr_bnds[model->getNLocalVars()];
         model->getObjGradient(obj);
         model->getColLowBounds(lwr_bnds);
         model->getColUpBounds(upr_bnds);
         double *op=obj;
         double *lb = lwr_bnds;
         double *ub = upr_bnds;
         for(int j=col_offset[model->getName()]; j<next_col; ++j,++op,++lb,++ub) {
            sparse_col &col = cols[j];
            if(*op!=0) col.add_entry(0, *op);
            if(*lb == *ub) { // Fixed
               bnds.push_back(mps_bound("FX", j, *lb));
            } else if(*lb==-inf && *ub==inf) {
               bnds.push_back(mps_bound("FR", j, 0.0));
            } else {
               if(*lb!=0.0)
                  bnds.push_back(mps_bound("LO", j, *lb));
               if(*ub!=inf)
                  bnds.push_back(mps_bound("UP", j, *ub));
            }
         }
      }

      // Name rows and identify constraint types
      {
         double lwr_bnds[model->getNLocalCons()];
         double upr_bnds[model->getNLocalCons()];
         model->getRowBounds(lwr_bnds, upr_bnds);
         double *lb = lwr_bnds;
         double *ub = upr_bnds;
         for(int j=row_offset[model->getName()]; j<next_row; ++j,++lb,++ub) {
            if(*lb==-inf && *ub==inf) {
               rows[j] = 'N';
            } else if(*lb==-inf) {
               rows[j] = 'L';
               if(*ub!=0) rhs[j] = *ub;
            } else if(*ub== inf) {
               rows[j] = 'G';
               if(*lb!=0) rhs[j] = *lb;
            } else if(*lb==*ub) {
               rows[j] = 'E';
               if(*lb!=0) rhs[j] = *lb;
            } else {
               rows[j] = 'L';
               if(*ub!=0) rhs[j] = *ub;
               ranges[j] = *ub-*lb;
            }
         }
      }
   }

   // Now aquire matrix information by ascending and descending tree
   // Recall that model->getJacobianOfIntersection(jm) returns us the
   // block in the current model's row using variables in block jm
   for(ExpandedModelInterface::child_iterator i=root->cbegin(); i!=root->cend(); ++i) {
      ExpandedModelInterface *model = *i;
      for(ExpandedModelInterface::child_iterator j=model->cbegin(); j!=model->cend(); ++j) {
         ExpandedModelInterface *jm = *j;
         if(model->getNzJacobianOfIntersection(jm)==0) continue;
         int nz = model->getNzJacobianOfIntersection(jm);
         int colbeg[jm->getNLocalVars()+1];
         int collen[jm->getNLocalVars()];
         int rownbs[nz];
         double elts[nz];
         model->getJacobianOfIntersection(jm, colbeg, collen, rownbs, elts);
         int offset = row_offset[model->getName()];
         for(unsigned int p=0; p<jm->getLocalVarNames().size(); ++p) {
            sparse_col &col = cols[col_offset[jm->getName()]+p];
            for(int k=colbeg[p]; k<colbeg[p]+collen[p]; ++k) {
               col.add_entry(rownbs[k]+offset, elts[k]);
            }
         }
      }

      for(ExpandedModelInterface::ancestor_iterator j=model->abegin(); j!=model->aend(); ++j) {
         ExpandedModelInterface *jm = *j;
         if(model->getNzJacobianOfIntersection(jm)==0) continue;
         int nz = model->getNzJacobianOfIntersection(jm);
         int colbeg[jm->getNLocalVars()+1];
         int collen[jm->getNLocalVars()];
         int rownbs[nz];
         double elts[nz];
         model->getJacobianOfIntersection(jm, colbeg, collen, rownbs, elts);
         int offset = row_offset[model->getName()];
         for(unsigned int p=0; p<jm->getLocalVarNames().size(); ++p) {
            sparse_col &col = cols[col_offset[jm->getName()]+p];
            for(int k=colbeg[p]; k<colbeg[p]+collen[p]; ++k) {
               col.add_entry(rownbs[k]+offset, elts[k]);
            }
         }
      }
   }

   string newfilename = "../" + filename;
   ofstream fout(newfilename.c_str());
   writeMps(fout, "Test", rows, rhs, ranges, bnds, cols);
   fout.close();
}

string mps_float(double d) {
   ostringstream oss;
   oss << setw(11) << d;
   if(oss.str().find('e') != oss.str().npos) {
      ostringstream oss2;
      if(d<1) { // Very small
         oss2 << fixed << setprecision(10) << d;
      } else { // Very big
         oss2 << fixed << setprecision(0) << d;
      }
      return oss2.str();
   }
   if(oss.str().find('.') == oss.str().npos) oss << '.';
   return oss.str();
}

string concat(char c, int n) {
   ostringstream oss;
   oss << c << n;
   return oss.str();
}

void writeMps(ostream &out, string name, map<int,char> rows,
      map<int,double> rhs, map<int,double> ranges, list<mps_bound> bnds,
      map<int,sparse_col> cols) {
   out << "NAME          " << name << endl;
   
   out << "ROWS" << endl;
   for(map<int,char>::iterator i=rows.begin(); i!=rows.end(); ++i) {
      out << " " << i->second << "  " << concat('r',i->first) << endl;
   }

   out << "COLUMNS" << endl;
   for(map<int,sparse_col>::iterator i=cols.begin(); i!=cols.end(); ++i) {
      sparse_col &col = i->second;
      bool rfirst = true;
      for(int j=0; j<col.size(); ++j) {
         if(rfirst) out << "    " << setw(8) << concat('c',i->first) << "  ";
         else out << "   ";
         out << setw(8) << concat('r', col.row[j]) << "  ";
         out << setw(12) << mps_float(col.val[j]);
         if(!rfirst) out << endl;
         rfirst = !rfirst;
      }
      if(!rfirst) out << endl;
   }

   out << "RHS" << endl;
   for(map<int,double>::iterator i=rhs.begin(); i!=rhs.end(); ++i) {
      out << "    rhs1      " << setw(8) << concat('r',i->first) << "  ";
      out << setw(12) << mps_float(i->second) << endl;
   }

   out << "RANGES" << endl;
   for(map<int,double>::iterator i=ranges.begin(); i!=ranges.end(); ++i) {
      out << "    range1    " << setw(8) << concat('r',i->first) << "  ";
      out << setw(12) << mps_float(i->second) << endl;
   }

   out << "BOUNDS" << endl;
   for(list<mps_bound>::iterator i=bnds.begin(); i!=bnds.end(); ++i) {
      out << " " << setw(2) << i->type << " ";
      out << setw(8) << "bnd1" << "  ";
      out << setw(8) << concat('c',i->id) << "  ";
      out << setw(12) << mps_float(i->value) << endl;
   }

   out << "ENDATA" << endl;
}
