$if not set jacFileName $abort --jacFileName missing
$GDXIN '%jacFileName%'
$onempty

Set i(*) Equation names ;
$LOADDC i<i.dim1

Set j(*) Variable names ;
$LOADDC j<j.dim1

Set jobj(j) Objective name ;
$LOADDC jobj

Scalar objcoef Objective coefficient ;
$LOADDC objcoef

Equation e(i) Equations ;
$LOADDC e

free     Variable x(j) Variables ;
$LOADDC x

Parameter A(i,j) Jacobian ;
$LOADDC A

Set ANl(i,j) Non-linear Jacobian indicator / /;

$setnames '%jacFileName%' fp fn fe
$log '%fp%%fn%_novenames%fe% 
$gdxout '%fp%%fn%_novenames%fe%'
$unload
