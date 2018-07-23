Set i equality-rows    / i1*i7 /
    k inequali-rows    / k1*k7 /
    j columns / j1*j9 /;

parameter g(j) obj coefficients / j1 2, j2 2, j3 2, j4 2, j5 2, j6 2, j7 2, j8 2, j9 2 /
          b(i) right hand side  / i1 3, i2 10, i3 5, i4 -12, i5 7, i6 4, i7 6 /
          cupp(k) right hand side  / k1 3, k2 22, k3 -2, k4 -1, k5 -0.5, k6 8, k7 10 /;

Table A(i,j)
    j1   j2    j3      j4    j5    j6    j7    j8    j9
i1   2    1
i2   3    7
*i3   2    3                   1
*i4  -4   -2                  -2
i3   2    1     1             1
i4  -4   -2            -2    -2
i5   6    1                 
i6   1                                    3
i7   1    1                         1     1     1     1
;

Table C(k,j)
    j1   j2   j3     j4    j5    j6    j7    j8    j9
k1   3
k2   6   13
k3  -4   -4                -1
k4  -2   -2                 1         
k5  -3   -1          
k6        4                       1     1     2
k7   1    1                             1     1
;

Positive Variables x(j)  / j1.lo 1, j3.up 10, j4.up 5, j5.up 40 /;
*makes it infeasible: j2.up 0.5
*j4.lo 1

Variable           z      objective variable
Equations          e(i)   equality equations
                   ie(k)  inequality equations
                   defobj objective function;

defobj.. z =e= sum(j, g(j)*x(j));
e(i)..   sum(j, A(i,j)*x(j)) =e= b(i);
ie(k)..  sum(j, C(k,j)*x(j)) =l= cupp(k);

Model m /all/ ;

$ifthen %METHOD%==PIPS
*annotations for variables:
  x.stage('j5') = 2;
  x.stage('j6') = 3;
  x.stage('j7') = 3;
  x.stage('j8') = 3;
  x.stage('j9') = 3;
*annotations for equations:
  e.stage('i1') = 1;
  e.stage('i2') = 1;
  e.stage('i3') = 2;
  e.stage('i4') = 2;
  e.stage('i5') = 2;
  e.stage('i6') = 3;
  e.stage('i7') = 4;
  ie.stage('k1') = 1;
  ie.stage('k2') = 1;
  ie.stage('k3') = 2;
  ie.stage('k4') = 2;
  ie.stage('k5') = 2;
  ie.stage('k6') = 3;
  ie.stage('k7') = 4;
  defobj.stage  = 4;


*For creation of gdx files:
$ echo jacobian exampleAC_parallelRow.gdx > convertd.opt
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


