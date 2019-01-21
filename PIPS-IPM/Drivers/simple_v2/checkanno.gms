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

$if not set maxStages $set maxStages 1000
Set s stages /0*%maxStages%/;
Scalar maxVarStage; maxVarStage = smax(j, x.stage(j));
if (maxVarStage+1>%maxStages%,
  put_utility 'log' / 'Increase --maxStages to at least ' (maxVarStage+2):0:0 ' for a proper report';
);

parameter repSize(s,*) / /; option repSize:0:1:1;
loop(s$sameas('0',s), loop(i, repSize(s+e.stage(i),'m orig') = repSize(s+e.stage(i),'m orig')+1));
loop(s$sameas('0',s), loop(j, repSize(s+x.stage(j),'n orig') = repSize(s+x.stage(j),'n orig')+1));
$ifthen set reportOnly
  display repSize;
$ exit
$endif

* Remove first stage variables from Aij
Set Aji(j,i), Aij(i,j), j0(j);
option Aji<A;
j0(j) = x.stage(j)=1;
Aji(j0,i) = no;
option Aij<Aji;
option clear=Aji;

Parameter maxStage(i), minStage(i),  badAnnotation /0/;
maxStage(i) = smax(Aij(i,j), x.stage(j));
minStage(i) = smin(Aij(i,j), x.stage(j));

$ifthen not set skipCheck
  loop(i,
    if (mapval(maxStage(i))=mapval(-inf),
      if (e.stage(i)<>1,
        put_utility 'log' / i.tl:0 ' has stage ' e.stage(i):0:0 ' but should have stage 1';
        badAnnotation = badAnnotation+1;
      )
    elseif maxStage(i)=minStage(i),
      if (minStage(i)<>e.stage(i),
        put_utility 'log' / i.tl:0 ' has stage ' e.stage(i):0:0 ' but should have stage ' minStage(i):0:0;
        badAnnotation = badAnnotation+1;
      )
    elseif maxStage(i)>minStage(i),
      if (maxVarStage+1<>e.stage(i),
        put_utility 'log' / i.tl:0 ' has stage ' e.stage(i):0:0 ' but should have stage ' (maxVarStage+1):0:0;
        badAnnotation = badAnnotation+1;
      );      
    else
      put_utility 'log' / i.tl:0 ' unexpected min/max stage: ' minStage(i):0:0 ',' maxStage(i):0:0 ',' e.stage(i):0:0
      abort 'unexpected behavior';
    )
  );
  put_utility$badAnnotation 'log' / 'Original jacobian has no proper annotation';
$else
  badAnnotation = 1;
$endif

$ifthen not set skipFix
  if (badAnnotation,
    e.stage(i) = 1              $(mapval(maxStage(i))=mapval(-inf)) +
                 minStage(i)    $(maxStage(i)=minStage(i)) +
                 (maxVarStage+1)$(maxStage(i)>minStage(i));
$   setnames '%jacFileName%' fp fn fe
    put_utility 'log' / 'Writing modified jacobian with proper annotation';
    option gdxuels=full;
    execute_unload '%fp%%fn%_fixedanno_novenames%fe%', i, j, jobj, objcoef, e, x, A, ANl;
    loop(s$sameas('0',s), loop(i, repSize(s+e.stage(i),'m fixed') = repSize(s+e.stage(i),'m fixed')+1));
    loop(s$sameas('0',s), loop(j, repSize(s+x.stage(j),'n fixed') = repSize(s+x.stage(j),'n fixed')+1));
    display repSize;
  );
$endif

