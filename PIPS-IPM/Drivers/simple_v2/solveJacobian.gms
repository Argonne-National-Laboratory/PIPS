
Set
   i 'Equations index'
   j 'Variable index'
Singleton set
   jobj(j) 'Objective variable index';

Scalar objcoef 'Min or Max';

Equation e(i) 'Equation';

Variable x(j) 'Variables';

Parameter A(i,j) 'Jacobian'

$if not set jacfile $set jacfile allblocksPips.gdx
$gdxin %jacfile%
$load i j jobj objcoef e x A
$gdxIn

Set
   ie(i) '=e= equations' / #i /
   il(i) '=l= equations'
   ig(i) '=l= equations'
   in(i) '=n= equations'
;
il(i)    = mapval(e.lo(i))=7;
ig(i)    = mapval(e.up(i))=6;
ie(il)   = no;
ie(ig)   = no;
in(i)    = il(i) and ig(i);
il(in)   = no;
ig(in)   = no;

Parameter
   rhs(i)  'right hand side'
   AObj(i) 'objective coefficient in row i';

AObj(i)   = A(i,jobj);
A(i,jobj) = 0;
rhs(ie)   = e.lo(ie);
rhs(il)   = e.up(il);
rhs(ig)   = e.lo(ig);

* Model
Equation defie(i), defil(i), defig(i), defin(i);
Variable z;

defie(ie).. sum(j, x(j)*A(ie,j)) + z*AObj(ie) =e= rhs(ie);

defil(il).. sum(j, x(j)*A(il,j)) + z*AObj(il) =l= rhs(il);

defig(ig).. sum(j, x(j)*A(ig,j)) + z*AObj(ig) =g= rhs(ig);

defin(in).. sum(j, x(j)*A(in,j)) + z*AObj(in) =g= rhs(in);

Model fromJacobian / all - e /;

option limRow=0, limCol=0, solPrint=off, resLim=1e9;
fromJacobian.dictFile = 0;
fromJacobian.reslim = 60*60*24*5;

$iftheni %TARGET%==LP
  option lp=convert; fromJacobian.optfile=1;
$ echo CplexLP %JACFILE%.lp > convert.opt
$elseifi %TARGET%==MPS
  option lp = convert; fromJacobian.optfile=1;
$ echo CplexMPS %JACFILE%.mps > convert.opt
$elseifi %TARGET%==SOPLEXLP
  option lp=soplex, iterlim=0;
* integer3=2 triggers creation of lp file
  option integer3=2;
$elseifi %TARGET%==CPLEXLP
  option lp=cplex, iterlim=0; fromJacobian.optfile=1;
$ echo writelp %JACFILE%.lp > cplex.opt
$endif

if (objcoef<0,
   solve fromJacobian max z using lp;
else
   solve fromJacobian min z using lp;
)

