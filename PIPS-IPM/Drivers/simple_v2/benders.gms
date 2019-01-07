$if not set BENDERSMAXITER $ set BENDERSMAXITER    5000

Set kk                   benders cuts        / k1*k%BENDERSMAXITER% /
    k(kk)                active benders cuts
;
Parameters
    bconst(kk)           constant part in benders cut
    bcoef(kk,tt,rr1,rr2) flow coefficient in benders cut
    bcoef2(kk,rr,e)      emission_split coefficient in benders cut
;

Variable
    Z                    future cost
Equations
    bobj                 benders master problem objective
    bcut(kk)             benders cut
;

bobj..
    obj =e= sum(net, LINK_ADD_CAP(net) * cost_link_add(net)) + Z;

bcut(k)..
    Z =g= bconst(k) + sum((t,net), bcoef(k,t,net) * FLOW(t,net))
                    + sum((rr,e), bcoef2(k,rr,e) * EMISSION_SPLIT(rr,e));

model master / bobj, bcut, eq_emission_cap, eq_link_capacity  /;
model sub    / simple - eq_emission_cap - eq_link_capacity /;

master.optfile = 1;
sub.optfile    = 1;

* set CPLEX feasibility and optimality tolerance
$onecho > cplex.opt
names no
lpmethod 1
*eprhs 1e-9
*epopt 1e-9
$offecho

parameter
    rep                        report parameter
    lb                         lower bound                         / -inf /
    ub                         upper bound
    gub                        global upper bound                  / inf /
    done                       termination indicator               / 0 /
    best_flow(tt,rr1,rr2)      flow corresponding to gub           
    best_link_add_cap(rr1,rr2) link_add_cap corresponding to gub   
    best_emission_split(rr,e)  emission_split corresponding to gub / #rr.#e 0 /
;
best_flow(tt,net)      = 0;
best_link_add_cap(net) = 0;

option limrow=0, limcol=0, solprint=silent, solvelink=5;

* First iteration with y=0
t(tt)                   = yes;
FLOW.fx(t,net)          = 0;
LINK_ADD_CAP.fx(net)    = 0;
EMISSION_SPLIT.fx(rr,e) = sum(t,demand(t,rr)) / sum((rr2,t), demand(t,rr2));
r(rr)                   = no;

loop(kk$(not done and ord(kk)<%BENDERSMAXITER%),
  ub = sum(net, LINK_ADD_CAP.l(net) * cost_link_add(net));
  bcoef(kk,t,net) = 0;
  bcoef2(kk,rr,e) = 0;
  bconst(kk)      = 0;
  loop(rr,
    r(rr) = yes;
    solve sub min obj use lp;
    abort$(sub.modelstat > 1) 'sub not solved to optimality';
    ub = ub + obj.l;

    bconst(kk) = bconst(kk) + sum(t, eq_power_balance.m(t,rr)               * demand(t,rr))
                            + sum((t,rp(rr,p)), eq_plant_capacity.m(t,rp)   * plant_cap(t,rp)*avail(t,rp))
                            + sum(rp(rr,p), eq_total_plant_capacity.m(rp)   * total_plant_cap(rp))
                            + sum((t,rs(rr,s)), eq_storage_capacity.m(t,rs) * storage_cap(rs));

    bcoef(kk,t,net(rr,rr2)) = bcoef(kk,t,net) + eq_power_balance.m(t,rr);
    bcoef(kk,t,net(rr2,rr)) = bcoef(kk,t,net) - eq_power_balance.m(t,rr) * link_efficiency(t,net);
    bcoef2(kk,rr,e)         = bcoef2(kk,rr,e) + eq_emission_region.m(rr,e) * total_emission_cap(e);

    r(rr) = no;
  );
  if(ub<gub,
    best_flow(t,net)          = FLOW.l(t,net);
    best_link_add_cap(net)    = LINK_ADD_CAP.l(net);
    best_emission_split(rr,e) = EMISSION_SPLIT.l(rr,e);
    gub = ub;
  );

  FLOW.lo(t,net)          = 0;
  FLOW.up(t,net)          = inf;
  LINK_ADD_CAP.lo(net)    = 0;
  LINK_ADD_CAP.up(net)    = link_max_add_cap(net);
  EMISSION_SPLIT.lo(rr,e) = 0;
  EMISSION_SPLIT.up(rr,e) = inf;

  k(kk) = yes;

  solve master min obj us lp;

  abort$(master.modelstat > 1) 'master not solved to optimality';

  FLOW.fx(t,net)          = FLOW.l(t,net);
  LINK_ADD_CAP.fx(net)    = LINK_ADD_CAP.l(net);
  EMISSION_SPLIT.fx(rr,e) = EMISSION_SPLIT.l(rr,e);

  lb               = OBJ.l;
  done             = (gub - lb)/abs(gub) < 1e-8;
  rep(kk,'lb')     = lb;
  rep(kk,'gloub')  = gub;
  rep(kk,'locub')  = ub;
);
* If converged, reestablish the best solution
if (done,
  FLOW.fx(t,net)          = best_flow(t,net);
  LINK_ADD_CAP.fx(net)    = best_link_add_cap(net);
  EMISSION_SPLIT.fx(rr,e) = best_emission_split(rr,e);
  loop(rr,
    r(rr) = yes;
    solve sub min obj use lp;
    abort$(sub.modelstat > 1) 'sub not solved to optimality';
    r(rr) = no;
  );
)
