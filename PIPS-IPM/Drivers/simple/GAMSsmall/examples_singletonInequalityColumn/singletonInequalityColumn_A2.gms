* ./gamsexample.sh -NP=3 -BLOCKS=4 -GAMSFILE=./singletonInequalityColumn_A2 -PRESOLVE=true


* if you wanna use this example please hack in useLinkStructure = true; in sData.C
Set i rows    / i1*i15 /
* subset of linking rows
    j columns / j1*j14 /;

parameter g(j) obj coefficients / j1 0, j2 1, j3 1, j4 1, j5 1, j6 1, j7 1, j8 1, j9 1, j10 1, j11 1, j12 1, j13 1, j14 1 /
          bA(i) right hand side  / i1 1, i2 -2, i3 0, i4 1, i5 -2, i6 0, i7 1, i8 -2, i9 0, i10 1, i11 -2, i12 0, i13 -10, i14 -3 /
*          clow(i) c left hand side    / i1 0, i2 1, i3 -1, i4 -100, i5 -100, i6 -100, i7 -100, i8 -100, i9 -100, i10 -100, i11 -100, i12 -100, i13 -100, i14 -100, i15 -100 /
          cupp(i) c right hand side   / i1 100, i2 1, i3 -1, i4 100, i5 100, i6 100, i7 100, i8 -1, i9 100, i10 -1, i11 100, i12 0, i13 100, i14 100, i15 100 /


* in this example the singletonColumnPresolver should substitute j1 (free) and put it into the objective
Table A(i,j)
     j1    j2    j3    j4    j5    j6    j7    j8    j9   j10   j11   j12   j13   j14
i1                3     2    -1
i2                2    -2     4
i3               -1   0.5    -1
i4                                  3     2    -1
i5                                  2    -2     4
i6                                 -1   0.5    -1
i7                                                    3     2    -1
i8                                                    2    -2     4
i9                                                   -1   0.5    -1
i10                                                                     3     2    -1
i11                                                                     2    -2     4
i12                                                                    -1   0.5    -1
i13               1     1     1     1     1     1     1     1           1     1     1
i14               1     1                 1
;       
*      2   -15    1    -2    -2     1    -2    -2     1    -2    -2     1    -2    -2
* expected values for x full determined by Ax=b

* the singletonInequalityRow i10:
* -31 x1 + x9 - 10 x10 - 20 x11 <= -1
* in the final solution x9 = 1, x10 = -2, x11 = -2 thus we get
* -31 x1 + 1 - 10 + 40 <= -1 <=> -31 x1 <= -62 <=> x1 <= 2 -> making the objective coeff negative will lead to x1 = 2 
* we also have a fixing to do i8:
* 0.5 j2 + j9 -3 j10 <= -1 -> in the final solution
* 0.5 j2 + 1 + 6 <= -1 
* 0.5 j2 <= -8 <=> 8 < 0.5 j2 <=> -16 >= j2 -> fix j2 at it's known lower bound (since the signums meet as well as it's objective direction)
* expected objective is -2
Table C(i,j)
    j1    j2    j3    j4    j5    j6    j7    j8    j9   j10   j11   j12   j13   j14
i1               1
i2               1     1     1
i3               1           1
i4                                 1     1
i5                                       1     1
i6                                       1     1
i7                                                   1     1     1
i8       0.5                                         1    -3
i9                                                   1           1
i10  -31                                             1   -10   -20
i11                                                                    1     1     1
i12                                                                          1     
i13                                                                                1
i14              1     1           1
i15              1           1     1
;      

Variables          x(j) / j2.lo -20 /
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
  x.stage('j6') = 2;
  x.stage('j7') = 2;
  x.stage('j8') = 2;
  x.stage('j9') = 3;
  x.stage('j10') = 3;
  x.stage('j11') = 3;
  x.stage('j12') = 4;
  x.stage('j13') = 4;
  x.stage('j14') = 4;
*annotations for equations:
  e.stage('i1') = 1;
  e.stage('i2') = 1;
  e.stage('i3') = 1;
  e.stage('i4') = 2;
  e.stage('i5') = 2;
  e.stage('i6') = 2;
  e.stage('i7') = 3;
  e.stage('i8') = 3;
  e.stage('i9') = 3;
  e.stage('i10') = 4;
  e.stage('i11') = 4;
  e.stage('i12') = 4;
  e.stage('i13') = 5;
  e.stage('i14') = 5;
*  ge.stage('i1') = 1;
*  ge.stage('i2') = 1;
*  ge.stage('i3') = 1;
*  ge.stage('i4') = 2;
*  ge.stage('i5') = 2;
*  ge.stage('i6') = 2;
*  ge.stage('i7') = 3;
*  ge.stage('i8') = 3;
*  ge.stage('i9') = 3;
*  ge.stage('i10') = 4;
*  ge.stage('i11') = 4;
*  ge.stage('i12') = 4;
*  ge.stage('i13') = 5;
*  ge.stage('i14') = 5;
  le.stage('i1') = 1;
  le.stage('i2') = 1;
  le.stage('i3') = 1;
  le.stage('i4') = 2;
  le.stage('i5') = 2;
  le.stage('i6') = 2;
  le.stage('i7') = 3;
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
$ echo jacobian singletonInequalityColumn_A2.gdx > convertd.opt
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
