* ./gamsexample.sh -NP=3 -BLOCKS=4 -GAMSFILE=./singletonEqualityColumn_multiple_resulting_noLink -PRESOLVE=true


* if you wanna use this example please hack in useLinkStructure = true; in sData.C
Set i rows    / i1*i21 /
    ii ineqrows / ii1*ii19 /
* subset of linking rows
    j columns / j1*j16 /;

parameter g(j) obj coefficients / j1 -10, j2 -10, j3 0, j4 0, j5 0, j6 -10, j7 0, j8 -5, j9 0, j10 0, j11 0, j12 0, j13 0, j14 1, j15 1, j16 1 /
          bA(i) right hand side  / i1 0, i2 1, i3 -2, i4 0, i5 1, i6 -2, i7 0, i8 0, i9 0, i10 1, i11 0, i12 -2, i13 0, i14 1, i15 -2, i16 0, i17 -10, i18 -5, i19 -8, i20 -3, i21 4 /
*          clow(ii) c left hand side   / ii1 0, ii2 1, ii3 -1, ii4 -100, ii5 -100, ii6 -100, ii7 -100, ii8 -100, ii9 -100, ii10 -100, ii11 -100, ii12 -100, ii13 -100, ii14 -100, ii15 -100, ii16 -100, ii17 -100, ii18 -100, ii19 -100/
          cupp(ii) c right hand side   / ii1 100, ii2 100, ii3 100, ii4 100, ii5 100, ii6 100, ii7 100, ii8 100, ii9 100, ii10 100, ii11 100, ii12 100, ii13 100, ii14 100, ii15 100, ii16 100, ii17 100, ii18 100, ii19 100/


* in this example the singletonColumnPresolver should substitute j1 (free) and put it into the objective
Table A(i,j)
    j1    j2    j3    j4    j5    j6    j7    j8    j9   j10   j11   j12   j13   j14   j15   j16
i1   1         0.1   0.1   0.1
i2               3     2    -1
i3               2    -2     4
i4              -1   0.5    -1
i5                                       3           2    -1
i6                                       2          -2     4
i7                                 1        -0.5  0.05   0.1   
i8                                     0.1     1  0.05
i9                                      -1         0.5    -1
i10                                                              3     2    -1
i11        1                                                   0.1   0.1   0.1
i12                                                              2    -2     4
i13                                                             -1   0.5    -1
i14                                                                                3     2    -1
i15                                                                                2    -2     4
i16                                                                               -1   0.5    -1
i17              1     1    1            1           1     1     1     1           1     1     1
i18              1     1                             1           1     1           1           1
i19              1     1    1            1           1     1     1     1           1           1
i20                    1                             1           1
i21              1                       1                       1                 1
; 
*                1    -2    -2           1          -2    -2     1    -2    -2     1    -2    -2
* expected values for x full determined by Ax=b

Table C(ii,j)
     j1    j2    j3    j4    j5    j6    j7    j8    j9    j10   j11   j12   j13   j14   j15   j16
ii1               1
ii2               1     1
ii3               1           1
ii4               1     1  
ii5                           1
ii6                                       1
ii7                                                         1
ii8                                                   1     1
ii9                                                               1     1
ii10                                                                    1
ii11                                                                          1
ii12                                                                                            1
ii13                                                                                1     1     
ii14                                                                                            1
ii15              1     1                 1                                                     1
ii16              1           1           1           1                 1
ii17              1     1                             1                       1
ii18              1           1           1           1                 1     1                 1
ii19              1     1                             1                       1           1
;

Variables          x(j) / j4.lo -5 /
Variable           z      objective variable
Equations          e(i)   equality equations
*                   ge(ii)  greater than inequality equations
                   le(ii)  less than inequality equations
                   defobj objective function;

defobj.. z =e= sum(j, g(j)*x(j));
e(i)..   sum(j, A(i,j)*x(j)) =e= bA(i);
*ge(ii)..  sum(j, C(ii,j)*x(j)) =g= clow(ii);
le(ii)..  sum(j, C(ii,j)*x(j)) =l= cupp(ii);

Model m /all/ ;


$ifthen %METHOD%==PIPS
*annotations for variables:
  x.stage('j6') = 2;
  x.stage('j7') = 2;
  x.stage('j8') = 2;
  x.stage('j9') = 2;
  x.stage('j10') = 2;
  x.stage('j11') = 3;
  x.stage('j12') = 3;
  x.stage('j13') = 3;
  x.stage('j14') = 4;
  x.stage('j15') = 4;
  x.stage('j16') = 4;
*annotations for equations:
  e.stage('i1') = 1;
  e.stage('i2') = 1;
  e.stage('i3') = 1;
  e.stage('i4') = 1;
  e.stage('i5') = 2;
  e.stage('i6') = 2;
  e.stage('i7') = 2;
  e.stage('i8') = 2;
  e.stage('i9') = 2;
  e.stage('i10') = 3;
  e.stage('i11') = 3;
  e.stage('i12') = 3;
  e.stage('i13') = 3;
  e.stage('i14') = 4;
  e.stage('i15') = 4;
  e.stage('i16') = 4;
  e.stage('i17') = 5;
  e.stage('i18') = 5;
  e.stage('i19') = 5;
  e.stage('i20') = 5;
  e.stage('i21') = 5;
*  ge.stage('ii1') = 1;
*  ge.stage('ii2') = 1;
*  ge.stage('ii3') = 1;
*  ge.stage('ii4') = 2;
*  ge.stage('ii5') = 2;
*  ge.stage('ii6') = 2;
*  ge.stage('ii7') = 3;
*  ge.stage('ii8') = 3;
*  ge.stage('ii9') = 3;
*  ge.stage('ii10') = 4;
*  ge.stage('ii11') = 4;
*  ge.stage('ii12') = 4;
*  ge.stage('ii13') = 5;
*  ge.stage('ii14') = 5;
  le.stage('ii1') = 1;
  le.stage('ii2') = 1;
  le.stage('ii3') = 1;
  le.stage('ii4') = 1;
  le.stage('ii5') = 1;
  le.stage('ii6') = 2;
  le.stage('ii7') = 2;
  le.stage('ii8') = 2;
  le.stage('ii9') = 3;
  le.stage('ii10') = 3;
  le.stage('ii11') = 3;
  le.stage('ii12') = 4;
  le.stage('ii13') = 4;
  le.stage('ii14') = 4;
  le.stage('ii15') = 5;
  le.stage('ii16') = 5;
  le.stage('ii17') = 5;
  le.stage('ii18') = 5;
  le.stage('ii19') = 5;
  defobj.stage  = 5;


*For creation of gdx files:
$ echo jacobian singletonEqualityColumn_multiple_resulting_noLink.gdx > convertd.opt
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
