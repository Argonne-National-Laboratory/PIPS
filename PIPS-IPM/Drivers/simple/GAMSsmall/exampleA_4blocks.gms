Set i rows    / i1*i10 /
* subset of linking rows
    linkRows(i) / i9, i10 /
* subset of inequality rows
*    ineqRows(i) / i1, i3, i5, i7 /
    j columns / j1*j9 /;

parameter g(j) obj coefficients / j1 2, j2 2, j3 2, j4 2, j5 2, j6 2, j7 2, j8 2, j9 2 /
          b(i) right hand side  / i1 3, i2 7, i3 8, i4 4, i5 10, i6 17, i7 14, i8 10, i9 3, i10 5 /;
*          cupp(ineqRows) right hand side  / i1 5, i3 5, i5 3, i7 3 /;

Table A(i,j)
    j1   j2   j3     j4    j5    j6    j7    j8    j9
i1   1    2
i2   3    4
i3   2         3      3
i4   1    2    1
i5   4    5                 1     
i6   6    7                 2     2
i7   1    1                             3     4     5
i8   1                                  2     3     4
i9   1    1    1
i10                   1     1     1     1           1
;

Positive Variables x(j)
Variable           z      objective variable
Equations          e(i)   equality equations
*                   ie(ineqRows)  inequality equations
                   defobj objective function;

defobj.. z =e= sum(j, g(j)*x(j));
e(i)..   sum(j, A(i,j)*x(j)) =e= b(i);
*ie(ineqRows)..  sum(j, C(ineqRows,j)*x(j)) =l= cupp(ineqRows);

Model m /all/ ;

$ifthen %METHOD%==PIPS
*annotations for variables:
  x.stage('j3') = 2;
  x.stage('j4') = 2;
  x.stage('j5') = 3;
  x.stage('j6') = 3;
  x.stage('j7') = 4;
  x.stage('j8') = 4;
  x.stage('j9') = 4;
*annotations for equations:
  e.stage('i1') = 1;
  e.stage('i2') = 1;
  e.stage('i3') = 2;
  e.stage('i4') = 2;
  e.stage('i5') = 3;
  e.stage('i6') = 3;
  e.stage('i7') = 4;
  e.stage('i8') = 4;
  e.stage('i9') = 5;
  e.stage('i10') = 5;
  defobj.stage  = 5;


*For creation of gdx files:
$ echo jacobian exampleA_4blocks.gdx > convertd.opt
  option lp=convertd;
  m.optfile = 1;
  solve m use lp min z;
$else
  option lp=cplex;
$ onecho > cplex.opt
  lpmethod 4
  solutiontype 2
$ offecho
  m.optfile = 1;
  solve m use lp min z;
$endif





