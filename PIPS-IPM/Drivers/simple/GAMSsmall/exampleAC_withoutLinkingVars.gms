Set i rows    / i1*i7 /
    j columns / j1*j8 /;

parameter g(j) obj coefficients / j1 2, j2 2, j3 2, j4 2, j5 2, j6 2, j7 2, j8 2 /
          b(i) right hand side  / i1 3, i2 10, i3 4, i4 6, i5 5, i6 1, i7 8 /;

Table A(i,j)
    j1   j2    j3      j4    j5    j6    j7    j8
i1   2    1
i2   3    7
i3              1       2           1
i4              3       2     1
i5                                        4     1
i6                                        1
i7   1    1     1       1     1     1     1     1
;

Positive Variables x(j)  / j3.up 10, j7.up 5 /;

Variable           z      objective variable
Equations          e(i)   equality equations
                   defobj objective function;

defobj.. z =e= sum(j, g(j)*x(j));
e(i)..   sum(j, A(i,j)*x(j)) =e= b(i);

Model m /all/ ;

$ifthen %METHOD%==PIPS
*annotations for variables:
  x.stage('j1') = 1;
  x.stage('j2') = 1;
  x.stage('j3') = 2;
  x.stage('j4') = 2;
  x.stage('j5') = 2;
  x.stage('j6') = 2;
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
  defobj.stage  = 4;

*For creation of gdx files:
$ echo jacobian exampleAC_withoutLinkingVars.gdx > convertd.opt
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


