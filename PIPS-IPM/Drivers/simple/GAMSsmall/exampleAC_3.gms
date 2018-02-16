Set i rows    / i1*i8 /
* subset of linking rows
    linkRows(i) / i7, i8 /
* subset of inequality rows
    ineqRows(i) / i1, i3, i5, i7 /
    j columns / j1*j8 /;

parameter g(j) obj coefficients / j1 1, j2 3, j3 5, j4 1, j5 2, j6 4, j7 2, j8 2 /
          bA(i) right hand side  / i1 2, i2 7, i3 2, i4 3, i5 1, i6 12, i7 1, i8 2 /
          clow(ineqRows) c left hand side    / i1 1, i3 3, i5 6, i7 2 /
          cupp(ineqRows) c right hand side   / i1 2, i3 3, i5 11, i7 3 /;

Table A(i,j)
    j1   j2   j3     j4    j5    j6    j7    j8
i1   2
i2        7
i3   2         2
i4                    3
i5        1                 1
i6        4                       3
i7   1         1
i8                    1     1                 1
;

Table C(ineqRows,j)
    j1   j2   j3     j4    j5    j6    j7    j8
i1   2         
i3   3         4
i5   5                      6                 1
i7   1    1    1            
;

Positive Variables x(j)
Variable           z      objective variable
Equations          e(i)   equality equations
                   ge(ineqRows)  greater than inequality equations
                   le(ineqRows)  less than inequality equations
                   defobj objective function;

defobj.. z =e= sum(j, g(j)*x(j));
e(i)..   sum(j, A(i,j)*x(j)) =e= bA(i);
ge(ineqRows)..  sum(j, C(ineqRows,j)*x(j)) =g= clow(ineqRows);
le(ineqRows)..  sum(j, C(ineqRows,j)*x(j)) =l= cupp(ineqRows);

Model m /all/ ;

$ifthen %METHOD%==PIPS
*annotations for variables:
  x.stage('j3') = 2;
  x.stage('j4') = 2;
  x.stage('j5') = 3;
  x.stage('j6') = 3;
  x.stage('j7') = 3;
  x.stage('j8') = 3;
*annotations for equations:
  e.stage('i1') = 1;
  e.stage('i2') = 1;
  e.stage('i3') = 2;
  e.stage('i4') = 2;
  e.stage('i5') = 3;
  e.stage('i6') = 3;
  e.stage('i7') = 4;
  e.stage('i8') = 4;
  ge.stage('i1') = 1;
  ge.stage('i3') = 2;
  ge.stage('i5') = 3;
  ge.stage('i7') = 4;
  le.stage('i1') = 1;
  le.stage('i3') = 2;
  le.stage('i5') = 3;
  le.stage('i7') = 4;
  defobj.stage  = 4;


*For creation of gdx files:
$ echo jacobian exampleAC_3.gdx > convertd.opt
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





