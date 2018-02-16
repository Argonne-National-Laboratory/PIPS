Set i rows    / i1,i2,i3,i4,i5,i6,i7,ii1,ii2,ii3,ii4,ii5,ii6,ii7 /
   eqrows(i)  / i1,i2,i3,i4,i5,i6,i7 /
   ineqrows(i)  / ii1,ii2,ii3,ii4,ii5,ii6,ii7 /
* subset of linking rows
    linkRows(eqrows) / i6, i7 /
    j columns / j1,j2,j3,j4,j5,j6,j7,j8 /;

parameter c(j) obj coefficients / j1 2, j2 2, j3 2, j4 2, j5 2, j6 2, j7 2, j8 2 /
          b(i) right hand side  / i1 2, i2 7, i3 3, i4 7, i5 7, i6 6, i7 4, ii1 2, ii2 7, ii3 3, ii4 7, ii5 7, ii6 6, ii7 4 /;

Table A(i,j)
    j1      j2       j3      j4   j5     j6      j7    j8
i1           1                                
i2           7                                  
i3           2        0.5                         
i4           16       400                  
i5                                        1   
i6           1        1      100                        1
i7   1                1                          1      0.1
ii1          1                                  
ii2          7                                  
ii3          2        512                         
ii4          16       1                  
ii5          512                          1  
ii6          1        1      50           0.1           0.1
ii7  0.2              1                          1      70
;



Positive Variables x(j) / j6.up 6 /;
Variable           z      objective variable
Equations          e(eqrows)
                   defobj objective function
                   f(ineqrows);

defobj.. z =e= sum(j, c(j)*x(j));
e(eqrows)..   sum(j, A(eqrows,j)*x(j)) =e= b(eqrows);
f(ineqrows)..   sum(j, A(ineqrows,j)*x(j)) =g= b(ineqrows);


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
  f.stage('ii3') = 2;
  f.stage('ii4') = 2;
  f.stage('ii5') = 3;
  f.stage('ii6') = 4;
  f.stage('ii7') = 4;

  defobj.stage  = 4;


*For creation of gdx files:
* $ echo jacobian epsBig.gdx > convertd.opt
$ echo jacobian epsBigXXX.gdx > convertd.opt
  option lp=convertd;
  m.optfile = 1;
  solve m use lp min z;
* begin new stuff
*  option lp = soplex;
*die folgende Zeile ist so eine Art hidden option um von Soplex das lp file zu bekommen
*  option integer3 = 2;
*  solve m use lp min z;
*end new stuff
$else
  option lp=cplex;
$ onecho > cplex.opt
  lpmethod 4
  solutiontype 2
$ offecho
  m.optfile = 1;
  solve m use lp min z;
$endif





