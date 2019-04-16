Set i rows    / i1*i10 /
* subset of linking rows
    linkRows(i) / i9, i10 /
    j columns / j1*j9 /;

parameter g(j) obj coefficients / j1 1, j2 1, j3 1, j4 1, j5 1, j6 1, j7 1, j8 1, j9 1 /
          bA(i) right hand side  / i1 2, i2 7, i3 -13, i4 5, i5 8, i6 16, i7 3, i8 1, i9 33, i10 7 /
          clow(i) c left hand side    / i1 -100, i2 -100, i3 -100, i4 -100, i5 -100, i6 -100, i7 -100, i8 -100, i9 -2, i10 -2 /
          cupp(i) c right hand side   / i1 100, i2 100, i3 100, i4 100, i5 100, i6 100, i7 100, i8 100, i9 4, i10 7 /

Table A(i,j)
    j1   j2   j3     j4    j5    j6    j7    j8    j9
i1   2
i2        7
i3  -2        -1    -11
i4   1    1           3
i5        2                 1     1
i6        4                -2     2
i7   1                              1e-16     1
i8                                     -1     1
i9   1    1    1      2     3     4     1     2 
i10  1         2      3     1           3
;
*    1    1    0      1     0     6     1     2     1
* expected values for x full determined by Ax=b

Table C(i,j)
    j1   j2   j3     j4    j5    j6    j7    j8    j9
i1   2
i2        7
i3   2         2
i4   1                3
i5        1                 1
i6        4                       3
i7   1         1
i8                                            1
i9   1         5      1   -11           1  -0.5     1
i10  1         5      1    12                       5
;

Positive Variables x(j)
Variable           z      objective variable
Equations          e(i)   equality equations
                   ge(i)  greater than inequality equations
                   le(i)  less than inequality equations
                   defobj objective function;

defobj.. z =e= sum(j, g(j)*x(j));
e(i)..   sum(j, A(i,j)*x(j)) =e= bA(i);
ge(i)..  sum(j, C(i,j)*x(j)) =g= clow(i);
le(i)..  sum(j, C(i,j)*x(j)) =l= cupp(i);

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
  ge.stage('i1') = 1;
  ge.stage('i2') = 1;
  ge.stage('i3') = 2;
  ge.stage('i4') = 2;
  ge.stage('i5') = 3;
  ge.stage('i6') = 3;
  ge.stage('i7') = 4;
  ge.stage('i8') = 4;
  ge.stage('i9') = 5;
  ge.stage('i10') = 5;
  le.stage('i1') = 1;
  le.stage('i2') = 1;
  le.stage('i3') = 2;
  le.stage('i4') = 2;
  le.stage('i5') = 3;
  le.stage('i6') = 3;
  le.stage('i7') = 4;
  le.stage('i8') = 4;
  le.stage('i9') = 5;
  le.stage('i10') = 5;
  defobj.stage  = 5;


*For creation of gdx files:
$ echo jacobian example_break_singletonrows.gdx > convertd.opt
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





