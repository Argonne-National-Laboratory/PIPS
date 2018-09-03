Set i rows    / i1,i2,i3,i4,i5,i6,i7 /
* subset of linking rows
    linkRows(i) / i6, i7 /
    j columns / j1,j2,j3,j4,j5,j6,j7,j8 /;

parameter c(j) obj coefficients / j1 2, j2 2, j3 2, j4 2, j5 2, j6 2, j7 2, j8 2 /
          b(i) right hand side  / i1 2, i2 7, i3 3, i4 7, i5 7, i6 6, i7 4 /;

Table A(i,j)
    j1      j2       j3     j4    j5     j6      j7    j8
i1   2                                          
i2           7                                  
i3   2                1                         
i4           5               2                  
i5           4                        0.4e-10   
i6   1       1        1      1            1             1
i7   1                1                          1      1
;

Positive Variables x(j)
Variable           z      objective variable
Equations          e(i)
                   defobj objective function;

defobj.. z =e= sum(j, c(j)*x(j));
e(i)..   sum(j, A(i,j)*x(j)) =g= b(i);

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
  e.stage('i3') = 2;
  e.stage('i4') = 2;
  e.stage('i5') = 3;
  e.stage(linkRows) = 4;
  defobj.stage  = 4;


*For creation of gdx files:
$ echo jacobian eps.gdx > convertd.opt
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





