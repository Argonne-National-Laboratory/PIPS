$if not set iterLim $set iterLim 10
$if not set resLim  $set resLim 200
$if not set rGap    $set rGap 0.1

option solvelink = 5, limrow = 0, limcol = 0, solprint = silent;

** Lagrangian Master Problem (D(lambda) approximation)
parameter
    lambda1(e),
    lambda2(tt,rr1,rr2),
    lambda3(rr1,rr2)
$if not set ll $set ll 0
$ifthen %ll%==1
    lambda2a(tt,rr1,rr2),
    lambda3a(rr1,rr2)
$endif
;

lambda1(e)          = 0;
lambda2(tt,rr1,rr2) = 0;
lambda3(rr1,rr2)    = 0;
$ifthen %ll%==1
lambda2a(tt,rr1,rr2) = 0;
lambda3a(rr1,rr2)    = 0;
$endif

* if RH used for preprocessing eq marginals are used for initial lambda
$ifthen %lgrhpp%==1
lambda1(e)          = rh_lambda1(e)          ;
lambda2(tt,rr1,rr2) = rh_lambda2(tt,rr1,rr2) ;
lambda3(rr1,rr2)    = rh_lambda3(rr1,rr2)    ;
$endif

$ontext
* use eq marginals of opt sol
lambda1(e)          = -eq_emission_cap.m(e);
lambda2(tt,rr1,rr2) = -eq_same_flow.m(tt,rr1,rr2);
lambda3(rr1,rr2)    = -eq_same_add_cap.m(rr1,rr2);
*display lambda1, lambda2, lambda3;
$offtext


Variable lagObj;
Equation defLagObj;
defLagObj..
    lagObj =e= sum(r, ROBJ(r))
             + sum(netx(r,net(rr1,rr2)), LINK_ADD_CAP(netx) * cost_link_add(net))
             + sum((e,r), lambda1(e)*EMISSION_SPLIT(r,e))

             - sum((t,netx(r,r,rr2)), lambda2(t,r,rr2)*FLOW(t,netx))
             + sum((t,netx(r,rr1,r)), lambda2(t,rr1,r)*FLOW(t,netx))

$ifthen %ll%==1
             + sum((t,netx(r,r,rr2)), lambda2a(t,r,rr2)*FLOW(t,netx))
             - sum((t,netx(r,rr1,r)), lambda2a(t,rr1,r)*FLOW(t,netx))
$endif
             - sum(netx(r,r,rr2), lambda3(r,rr2)*LINK_ADD_CAP(netx))
             + sum(netx(r,rr1,r), lambda3(rr1,r)*LINK_ADD_CAP(netx))

$ifthen %ll%==1
             + sum(netx(r,r,rr2), lambda3a(r,rr2)*LINK_ADD_CAP(netx))
             - sum(netx(r,rr1,r), lambda3a(rr1,r)*LINK_ADD_CAP(netx))
$endif
;

equation defLagHeurObj;
defLagHeurObj..
    lagObj =e= sum(r, ROBJ(r));

model simpleLag / simple - eq_obj - eq_emission_cap - eq_same_flow - eq_same_add_cap, defLagObj /;
model simpleLagHeur / simpleLag - defLagObj + defLagHeurObj /;

Set sg / 1*%iterLim% /;
Set sgact(sg);
Parameter
    subgradlam1(sg, e)              "lambda that defines sg'th subgradient"
    subgradlam2(sg, tt, rr1, rr2)   "ditto"
    subgradlam3(sg, rr1, rr2)       "ditto"
    subgradcoef1(sg, e)             "subgradient coefficients"
    subgradcoef2(sg, tt, rr1, rr2)  "ditto"
    subgradcoef3(sg, rr1, rr2)      "ditto"
    subgradDval(sg)                 "value of D(lambda)"
    proximity1(e)                   "proximity point"
    proximity2(tt, rr1, rr2)        "ditto"
    proximity3(rr1, rr2)            "ditto"
    proximityDval                   "value of D(proximity)"
$ifthen %ll%==1
    subgradlam2a(sg, tt, rr1, rr2)   "ditto"
    subgradlam3a(sg, rr1, rr2)       "ditto"
    proximity2a(tt, rr1, rr2)        "ditto"
    proximity3a(rr1, rr2)            "ditto"
$endif
    u                               "weight of proximity term";
Positive Variable
    lam1(e)                         "lambda"
$ifthen not %ll%==1
Variable
    lam2(tt, rr1, rr2)              "ditto"
    lam3(rr1, rr2)                  "ditto"
$else
Positive Variable
    lam2(tt, rr1, rr2)              "ditto"
    lam3(rr1, rr2)                  "ditto"
    lam2a(tt, rr1, rr2)              "ditto"
    lam3a(rr1, rr2)                  "ditto"
$endif

Variable dapprox                    "approximated value of D(lambda)";
Equation dapproxdef(sg)             "definition of D(.) approximation";
Variable bundleobj                  "objective variable"
Equation bundleobjdef               "objective function of bundle model";

dapproxdef(sgact(sg))..
    dapprox =l=   subgradDval(sg)
                + sum(e,                subgradcoef1(sg,e)     * (lam1(e)     - subgradlam1(sg,e)))
                + sum((t,net(rr1,rr2)), subgradcoef2(sg,t,net) * (lam2(t,net) - subgradlam2(sg,t,net)))
                + sum(net(rr1,rr2),     subgradcoef3(sg,net)   * (lam3(net)   - subgradlam3(sg,net)))
$ifthen %ll%==1
                - sum((t,net(rr1,rr2)), subgradcoef2(sg,t,net) * (lam2a(t,net) - subgradlam2a(sg,t,net)))
                - sum(net(rr1,rr2),     subgradcoef3(sg,net)   * (lam3a(net)   - subgradlam3a(sg,net)))
$endif
;

bundleobjdef..
    bundleobj =e= dapprox - u*(   sum(e,                sqr(proximity1(e)     - lam1(e)))
                                + sum((t,net(rr1,rr2)), sqr(proximity2(t,net) - lam2(t,net)))
                                + sum(net(rr1,rr2),     sqr(proximity3(net)   - lam3(net)))
$ifthen %ll%==1
                                + sum((t,net(rr1,rr2)), sqr(proximity2a(t,net) - lam2a(t,net)))
                                + sum(net(rr1,rr2),     sqr(proximity3a(net)   - lam3a(net)))
$endif
                               );
proximity1(e)              = lambda1(e);
proximity2(t,net(rr1,rr2)) = lambda2(t,net);
proximity3(net(rr1,rr2))   = lambda3(net);
$ifthen %ll%==1
proximity2a(t,net(rr1,rr2)) = lambda2a(t,net);
proximity3a(net(rr1,rr2))   = lambda3a(net);
$endif
proximityDval = -inf;

model Dfunc / dapproxdef, bundleobjdef /;

scalar violcnt, simplecost, hcost, converged /0/, lb, maxlb /-inf/, minub /inf/, predict /inf/;
$ifthen %lgrhpp%==1
minub = rh_obj;
$endif

scalar alpha  'subgradient method stepsize';
file lgfx /''/;
alias(sg,sgp);

parameter rep, rtol;
loop(sgp$(not converged),
* evaluate D(lambda), with Lagrangian subproblems
  lb = 0;
  simplecost = 0;
  loop(rr,
    r(rr) = yes;
    solve simpleLag min lagObj using lp;

    if (simpleLag.modelstat<>1, abort 'problems solving subproblem', r);
    lb = lb + lagObj.l;
    // "/2" is very subtle, but we need some way to count the copied variables just once
    simplecost = simplecost + ROBJ.l(rr) + sum(netx(rr,net(rr1,rr2)), LINK_ADD_CAP.L(netx)/2 * cost_link_add(net));
    r(rr) = no;
  );
  lb = lb - sum(e, lambda1(e)*1); // subtract the split RHS of 1
  maxlb$(lb>maxlb) = lb;  display lb;
*  u$sameas('1',sgp) = abs(lb)/2000;
*  u$sameas('1',sgp) = 1e2;
*  u$(mod(ord(sgp),10)=0) = u/2;
  u$sameas('1',sgp) = 1e-1;
*  u = 0;
* add new subgradient to bundle
  subgradDval(sgp) = lb;
  subgradlam1(sgp, e)              = lambda1(e);
  subgradlam2(sgp, t,net(rr1,rr2)) = lambda2(t,net);
  subgradlam3(sgp, net(rr1,rr2))   = lambda3(net);
$ifthen %ll%==1
  subgradlam2a(sgp, t,net(rr1,rr2)) = lambda2a(t,net);
  subgradlam3a(sgp, net(rr1,rr2))   = lambda3a(net);
$endif
  subgradcoef1(sgp, e)              = sum(rr, EMISSION_SPLIT.l(rr,e)) - 1;

  subgradcoef2(sgp, t,net(rr1,rr2)) = - FLOW.l(t,rr1,rr1,rr2) + FLOW.l(t,rr2,rr1,rr2);
  subgradcoef3(sgp, net(rr1,rr2))   = - LINK_ADD_CAP.l(rr1,rr1,rr2) + LINK_ADD_CAP.l(rr2,rr1,rr2);
  sgact(sgp) = yes;

* Update proximity, if current lambda gives better objective (D(lambda)) than value of current proximity point (D(proximity))
  if( subgradDval(sgp) > proximityDval,
     proximity1(e)              = lambda1(e);
     proximity2(t,net(rr1,rr2)) = lambda2(t,net);
     proximity3(net(rr1,rr2))   = lambda3(net);
$ifthen %ll%==1
     proximity2a(t,net(rr1,rr2)) = lambda2a(t,net);
     proximity3a(net(rr1,rr2))   = lambda3a(net);
$endif
     proximityDval = subgradDval(sgp);
     putclose lgfx '#### updating proximity point, new D(proximity)=', proximityDval /;
  );
  violcnt =   sum(e$(sum(rr,        EMISSION_SPLIT.l(rr,e)) > 1),1)
            + sum((t,net(rr1,rr2))$(abs(FLOW.l(t,rr1,rr1,rr2) - FLOW.l(t,rr2,rr1,rr2)) > eps), 1)
            + sum(net(rr1,rr2)$(abs(LINK_ADD_CAP.l(rr1,rr1,rr2) - LINK_ADD_CAP.l(rr2,rr1,rr2)) > eps), 1);
  hcost = 0;
*$ONTEXT
  if (violcnt = 0,
    minub$(simplecost<minub) = simplecost;
* Heuristic
  else
    if (ord(sgp)=1 or mod(ord(sgp),50)=0,
      hcost = 0;
      simpleLag.modelstat = 1;
      lambda1(e)          = 0        ;
      lambda2(tt,rr1,rr2) = 0;
      lambda3(rr1,rr2)    = 0   ;
$ifthen %ll%==1
      lambda2a(tt,rr1,rr2) = 0;
      lambda3a(rr1,rr2)    = 0  ;
$endif

      EMISSION_SPLIT.fx(rr,e)           = EMISSION_SPLIT.l(rr,e) / max(1,sum(rr2, EMISSION_SPLIT.l(rr2,e))) + max(0, (1-sum(rr2, EMISSION_SPLIT.l(rr2,e)))/card(rr2));
      FLOW.fx(t,netx(rr,rr1,rr2))       = min( FLOW.l(t,rr,rr1,rr2) , FLOW.l(t,rr,rr2,rr1) );
      LINK_ADD_CAP.fx(netx(rr,rr1,rr2)) = min( LINK_ADD_CAP.l(rr,rr1,rr2) , LINK_ADD_CAP.l(rr,rr2,rr1) );

      loop(rr$(simpleLag.modelstat=1),
        r(rr) = yes;
        solve simpleLag min lagObj using lp;

        hcost = hcost + ROBJ.l(rr) + sum(netx(rr,net(rr1,rr2)), LINK_ADD_CAP.L(netx)/2 * cost_link_add(net));
        rep(sgp,rr, 'robj') = ROBJ.l(rr);
        rep(sgp,rr, 'cap') = sum(netx(rr,net(rr1,rr2)), LINK_ADD_CAP.L(netx)/2 * cost_link_add(net))  ;
        r(rr) = no;
      );
      option clear = EMISSION_SPLIT, clear =  FLOW;
      LINK_ADD_CAP.lo(netx(rr,net)) = 0 ;
      LINK_ADD_CAP.up(netx(rr,net)) = link_max_add_cap(net)  ;
      hcost$(simpleLag.modelstat<>1) = inf;
      minub$(hcost<minub) = hcost;
    );
  );
*$OFFTEXT

* Update lambda
  dapprox.up = minub;
  option qcp=cplex;
  solve Dfunc max bundleobj using QCP;

  lambda1(e)          = lam1.l(e)         ;
  lambda2(tt,rr1,rr2) = lam2.l(tt,rr1,rr2);
  lambda3(rr1,rr2)    = lam3.l(rr1,rr2)   ;
$ifthen %ll%==1
  lambda2a(tt,rr1,rr2) = lam2a.l(tt,rr1,rr2);
  lambda3a(rr1,rr2)    = lam3a.l(rr1,rr2)   ;
$endif
  rtol = (minub-maxlb)/(abs(minub)+1e-3);
  putclose lgfx '##### ' sgp.tl:3:0 '    lb=', lb, ' predict=', predict, ' maxlb=', maxlb, ' minub=', minub, ' simplecost=', simplecost, ' violcnt=', violcnt, ' hcost=', hcost, ' gap=', rtol:10:4 /;
*abort 'sjkc';
  predict = dapprox.l;
  converged = (rtol<%rgap%) or (timeelapsed>%resLim%);
);
