* ./gamsexample.sh -NP=3 -BLOCKS=4 -GAMSFILE=nearlyParallelEqualityRowsBothSingletons_B0A2 -PRESOLVE=true

Set i rows    / i1*i18 /
* subset of linking rows
    j columns / j1*j18 /;

parameter g(j) obj coefficients / j1 0, j2 1, j3 1, j4 -1, j5 1, j6 1, j7 1, j8 1, j9 1, j10 1, j11 1, j12 1, j13 1 /
          bA(i) right hand side  / i1 -2, i2 1, i3 0, i4 -2, i5 0, i6 0, i7 -2, i8 0, i9 5, i10 -10, i11 -5,
            i12 -59, i13 -3, i14 1, i15 1, i16 -4, i17 -4, i18 -2 /
*          clow(i) c left hand side    / i1 0, i2 1, i3 -1, i4 -100, i5 -100, i6 -100, i7 -100, i8 -100, i9 -100, i10 -100, i11 -100, i12 -100, i13 -100, i14 -100, i15 -100 /
          cupp(i) c right hand side   / i1 100, i2 100, i3 100, i4 100, i5 100, i6 100, i7 100, i8 100, i9 100, i10 100, i11 100, i12 100, i13 100, i14 100, i15 100 /


* in this example the parallel row presolver should detect three parallel equality row pairs - one in B0 and two in B2 and delete them from the system
* for the one in B0 bounds get propagated but both variables are actually zero in the final solution
* for the one in B2 the propagated bound is tight
Table A(i,j)
    j1    j2    j3    j4    j5    j6    j7    j8    j9   j10   j11   j12   j13   j14   j15   j16   j17   j18
i1         1                      -6    -4     2
i2  -1                             3     2    -1
i3
i4                                 2    -2     4
i5                                -1   0.5    -1
i6                                 1     1           3     2    -1
i7                                                   2    -2     4
i8                                                  -1   0.5    -1
i9              10                                                           3     2    -1
i10                                                                   -2   -30   -20    10
i11                    1                 1                                  -1   0.5    -1
i12                         -2          10                                 -10     5   -10
i13                                1           1                             2    -2     4
i14                                                                                            3     2    -1
i15                                1    -1                                                     2    -2     4
i16                                      1     1                                              -1   0.5    -1
i17                                                  1     1     1           1     1
i18                                1     1                 1                                   1
;
*    0     0   0.4    -3   19.5    1    -2    -2     1    -2    -2     0     1    -2    -2     1    -2    -2
* expected values for x full determined by Ax=b

Table C(i,j)
    j1    j2    j3    j4    j5    j6    j7    j8    j9   j10   j11   j12   j13   j14   j15   j16   j17   j18
i1
i2                                 1
i3                                       1  
i4                                             1
i5                                                   1
i6                                                         1     
i7                                                               1
i8                                                                           1     1     1
i9                                                                           1
i10                                                                                1
i11                                                                                                  1     1
i12                                                                                            1     
i13                                                                                            1     1     1
i14                                1     1                 1                       1
i15                                1                       1                 1
;

Variables          x(j) / j1.lo -10, j1.up 5, j2.lo -1, j2.up 1, j3.lo -1, j12.lo 0, j12.up 30, j5.lo 19.5, j5.up 22 /
Variable           z      objective variable
Equations          e(i)   equality equations
*                   ge(i)  greater than inequality equations
                   le(i)  less than inequality equations
                   defobj objective function;

defobj.. z =e= sum(j, g(j)*x(j));
e(i)..   sum(j, A(i,j)*x(j)) =e= bA(i);
*ge(i)..  sum(j, C(i,j)*x(j)) =g= clow(i);
le(i)..  sum(j, C(i,j)*x(j)) =l= cupp(i);

Model m /all/ ;


$ifthen %METHOD%==PIPS
*annotations for variables:
  x.stage('j9') = 2;
  x.stage('j10') = 2;
  x.stage('j11') = 2;
  x.stage('j12') = 3;
  x.stage('j13') = 3;
  x.stage('j14') = 3;
  x.stage('j15') = 3;
  x.stage('j16') = 4;
  x.stage('j17') = 4;
  x.stage('j18') = 4;
*annotations for equations:
  e.stage('i1') = 1;
  e.stage('i2') = 1;
  e.stage('i3') = 1;
  e.stage('i4') = 1;
  e.stage('i5') = 1;
  e.stage('i6') = 2;
  e.stage('i7') = 2;
  e.stage('i8') = 2;
  e.stage('i9') = 3;
  e.stage('i10') = 3;
  e.stage('i11') = 3;
  e.stage('i12') = 3;
  e.stage('i13') = 3;
  e.stage('i14') = 4;
  e.stage('i15') = 4;
  e.stage('i16') = 4;
  e.stage('i17') = 5;
  e.stage('i18') = 5;
  le.stage('i1') = 1;
  le.stage('i2') = 1;
  le.stage('i3') = 1;
  le.stage('i4') = 1;
  le.stage('i5') = 2;
  le.stage('i6') = 2;
  le.stage('i7') = 2;
  le.stage('i8') = 3;
  le.stage('i9') = 3;
  le.stage('i10') = 3;
  le.stage('i11') = 4;
  le.stage('i12') = 4;
  le.stage('i13') = 4;
  le.stage('i14') = 5;
  le.stage('i15') = 5;
  defobj.stage  = 5;


*For creation of gdx files:
$ echo jacobian nearlyParallelEqualityRowsBothSingletons_B0A2.gdx > convertd.opt
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
