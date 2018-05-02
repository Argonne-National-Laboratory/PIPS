Set i rows    / i1*i6 /
    k ineqalityRows / k1*k9 /
    l geRows / k1 /
    j columns / j1*j8 /;

parameter g(j) obj coefficients / j1 2, j2 2, j3 2, j4 2, j5 2, j6 2, j7 2, j8 2 /
          b(i) right hand side  / i1 3, i2 10, i3 3, i4 6, i5 4, i6 9 /
          cupp(k) right hand side  / k1 3, k2 3, k3 6, k4 -1, k5 5, k6 7, k7 9, k8 10, k9 11/;

Table A(i,j)
    j1   j2    j3      j4    j5    j6    j7    j8
i1   2    1
i2   3    7
i3        1     1       1 
i4   5    1              
i5   1                              3
i6   1    1     2       1     1     1     1     1
;

Table C(k,j)
    j1   j2   j3     j4    j5    j6    j7    j8
k1   3
k2   3
k3   2         4      
k4  -2        -4
k5   1                            2
k6        4                       1     2
k7        4                       1     2
k8   1    1    1                  1     1
k9   1    1    1                  1     1
;

Table D(l,j)
    j1   j2   j3     j4    j5    j6    j7    j8
k1   3
;

Positive Variables x(j)  / j1.lo 1, j3.up 10, j7.up 5 /;
*makes it infeasible: j2.up 0.5

Variable           z      objective variable
Equations          e(i)   equality equations
                   ie(k)  inequality equations
                   ge(l)  greater than inequalities
                   defobj objective function;

defobj.. z =e= sum(j, g(j)*x(j));
e(i)..   sum(j, A(i,j)*x(j)) =e= b(i);
ie(k)..  sum(j, C(k,j)*x(j)) =l= cupp(k);
ge(l).. sum(j, D(l,j)*x(j)) =g= 1;

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
  e.stage('i6') = 4;
  ie.stage('k1') = 1;
  ie.stage('k2') = 1;
  ie.stage('k3') = 2;
  ie.stage('k4') = 2;
  ie.stage('k5') = 3;
  ie.stage('k6') = 3;
  ie.stage('k7') = 3;
  ie.stage('k8') = 4;
  ie.stage('k9') = 4;
  ge.stage('k1') = 1;
  defobj.stage  = 4;


*For creation of gdx files:
$ echo jacobian exampleAC_duplicateRow.gdx > convertd.opt
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


