$ eolcom //

$if  not set FROM                 $set FROM              0
$if  not set TO                   $set TO                1
$if  not set RESOLUTION           $set RESOLUTION        1
$if  not set NBREGIONS            $set NBREGIONS         4
$if  not set METHOD               $set METHOD            standard_lp
$if  not set SCALING              $set SCALING           0
$if  not set NOSLACK              $set NOSLACK           0

$ife %FROM%>%TO%         $abort 'FROM > TO'
$ife %RESOLUTION%<0      $abort 'Negative RESOLUTION forbidden'

$eval NBTIMESTEPS ceil((%TO%-%FROM%)*365*24/%RESOLUTION%)
$setnames "%gams.input%" SIMPLEDIR fn fe

* if there is no input data for the current settings available call data generator
$  if not exist  "%SIMPLEDIR%simple_data_%NBREGIONS%.gdx" $call gams  "%SIMPLEDIR%simple_data_gen.gms" --NBREGIONS=%NBREGIONS% lo=%gams.lo%
$  if errorlevel 1 $abort 'problems with automated data generation'



*   basic sets
Set rr                          'regions                                       '
    rrUel
    p                           'plants                                        '
    s                           'storages                                      '
    tt                          'time steps                                    ' / t1*t%NBTIMESTEPS% /
    ttUel
    ttX                         'time steps in input data                      '
    e                           'emissions                                     '
    type                        'plant type                                    '
    bb                          'blocks                                        '
;

Alias(rr,rr1,rr2,rr3), (tt,tt1,tt2);

*   Subsets and mappings
set r(rr)                       'active region                                 '
    t(tt)                       'subset of active time steps                   '
    t_fix(tt)                   'subset of time steps to be fixed after        '
    pr(p,rr)                    'plant to region mapping                       '
    rp(rr,p)                    'region to plant mapping                       '
    re(type)                    'renewable energy types                        '
    ptype(rr,p,type)            'plant type mapping                            '
    sr(s,rr)                    'storage to region mapping                     '
    rs(rr,s)                    'region to storage mapping                     '
    net(rr1,rr2)                'transmission links                            '
    netx(rr,rr1,rr2)            'copy of transmission links                    '
    time_map(tt,ttX)            'mapping of model to input time steps          '
    tX(ttX)                     'dynamic subset of data time steps             '
    tlast(tt)                   'last time step                                '
    tXlast(ttX)                 'last hour in base data                        '
;
Alias(r,r1,r2),(t,t1,t2);

Parameter
    plant_capX(rr,p,ttX)         'plant capacity                           [GW]'
    cost_unserved_demandX(ttX)   'price for unserved demand          [MEUR/GWh]'
    link_capX(rr1,rr2,ttX)       'transmission link capacity per hour     [GWh]'
    link_efficiencyX(rr1,rr2,ttX)'transmission link efficiency factor          '
    demandX(rr,ttX)              'demand for region per time step         [GWh]'
*   time parameters
    start(tt)                    'start hour of model time step                '
    end(tt)                      'end hour of model time step                  '
    startX(ttX)                  'start hour of data time step                 '
    endX(ttX)                    'end hour of data time step                   '
    overlap(tt,ttX)              'overlap of model and input time step         '
*   model parameters
    plant_cap(tt,rr,p)           'plant capacity in time step t           [GWh]'
    plant_cap2(rr,p,tt)          'plant capacity in time step t           [GWh]'
    yearly_plant_cap(rr,p)       'yearly plant capacity                   [GWh]'
    total_plant_cap(rr,p)        'plant capacity over total time span     [GWh]'
    plant_emission(rr,p,e)       'plant emission                     [tons/GWh]'
    cost_power_generation(rr,p)  'electricity production cost        [MEUR/GWh]'
    yearly_emission_cap(e)       'yearly emission cap                    [tons]'
    total_emission_cap(e)        'emission cap for total time span       [tons]'
    cost_emission(e)             'emission costs                     [MEUR/ton]'
    storage_cap(rr,s)            'storage capacity                        [GWh]'
    storage_efficiency(rr,s)     'effieciency factor of storage                '
    storage_efficiency_in(rr,s)  'effieciency factor of storage inflow         '
    storage_efficiency_out(rr,s) 'effieciency factor of storage outflow        '
    storage_max_in(rr,s)         'maximum inflow into storage             [GWh]'
    storage_max_out(rr,s)        'maximum outflow into storage            [GWh]'
    cost_unserved_demand(tt)     'price for unserved demand          [MEUR/GWh]'
    link_cap(tt,rr1,rr2)         'transmission link capacity per time step[GWh]'
    link_cap2(rr1,rr2,tt)        'transmission link capacity per time step[GWh]'
    link_efficiency(tt,rr1,rr2)  'transmission link efficiency factor          '
    link_efficiency2(rr1,rr2,tt) 'transmission link efficiency factor          '
    demand(tt,rr)                'demand for region per time step         [GWh]'
    demand2(rr,tt)               'demand for region per time step         [GWh]'
    plant_max_add_cap(rr,p)      'max additional plant capacity            [GW]'
    storage_max_add_cap(rr,s)    'max additional storage capacity         [GWh]'
    link_max_add_cap(rr1,rr2)    'max additional arc capacity             [GWh]'
    cost_plant_add(rr,p)         'cost for additional plant capacity  [MEUR/GW]'
    cost_storage_add(rr,s)       'cost for additional storage cap    [MEUR/GWh]'
    cost_link_add(rr1,rr2)       'cost for additional link capacity  [MEUR/GWh]'
    type_mult(type)              'generation cost multiplier                   '
    availX(ttX,rr,p)
    avail(tt,rr,p)
    avail2(rr,p,tt)
;

scalar
    tmpStart                     'helper storing time step start times         '
    tmpEnd                       'helper storing time step end times           '
;

* read data from gdx
$gdxin  "%SIMPLEDIR%simple_data_%NBREGIONS%.gdx"
$load ttX
$load rr
$load p s e pr rp type ptype re sr rs net
$load plant_capX cost_unserved_demandX link_capX link_efficiencyX demandX
$load yearly_plant_cap plant_emission cost_power_generation yearly_emission_cap
$load cost_emission storage_cap storage_efficiency storage_efficiency_in
$load storage_efficiency_out storage_max_in storage_max_out
$load plant_max_add_cap storage_max_add_cap link_max_add_cap cost_plant_add
$load cost_storage_add cost_link_add availX
$gdxin

tlast(tt)  = tt.last;
* compute mapping and overlap of model time steps to input time steps
start(tt)   = %FROM%*365*24 + (ord(tt)-1) * %RESOLUTION%;
end(tt)     = %FROM%*365*24 +  ord(tt)    * %RESOLUTION%;
end(tlast)  = %TO%*365*24;
startX(ttX) = ord(ttX)-1;
endX(ttX)   = ord(ttX);
tX(ttX)     = no;

loop(ttX$(not card(tX)), tX(ttX) = endX(ttX) > %FROM%*365*24;);
tmpStart = sum(tX, startX(tX));
tmpEnd   = sum(tX, endX(tX));

loop(tt,
  time_map(tt,tX) = yes;
  While(end(tt) > tmpEnd and tmpEnd < %TO%*365*24,
    tX(ttX)  = tX(ttX-1);
    tmpStart = sum(tX, startX(tX));
    tmpEnd   = sum(tX, endX(tX));
    time_map(tt,tX) = yes;
  );
  if(end(tt) = tmpEnd and not tlast(tt),
    tX(ttX)  = tX(ttX-1);
    tmpStart = sum(tX, startX(tX));
    tmpEnd   = sum(tX, endX(tX));
  );
  // if end(t) < tmpEnd we do not need to do anything
);
tXlast(tX) = tX.last;
abort$(not sum(time_map(tlast,tXlast),1)) 'Last time steps of base and model data do not map.';
overlap(time_map(tt,ttX)) = min(1, abs(min(0, max(end(tt),endX(ttX)) - min(start(tt),startX(ttX)) - (%RESOLUTION%+1))));

overlap(time_map(tlast,ttX)) = min(end(tlast),endX(ttX)) - max(start(tlast),startX(ttX));


* Compute parameters according to time span and resolution
plant_cap2(rp(rr,p),tt)      = sum(time_map(tt,ttX), plant_capX(rp,ttX) * overlap(tt,ttX));
total_plant_cap(rp(rr,p))    = (%TO% - %FROM%) * yearly_plant_cap(rp);
total_emission_cap(e)        = (%TO% - %FROM%) * yearly_emission_cap(e);
storage_efficiency(rs(rr,s)) = rPower(storage_efficiency(rs),%RESOLUTION%);
storage_max_in(rs(rr,s))     = %RESOLUTION% * storage_max_in(rs);
storage_max_out(rs(rr,s))    = %RESOLUTION% * storage_max_out(rs);
cost_unserved_demand(tt)     = smax(time_map(tt,ttX), cost_unserved_demandX(ttX));
link_cap2(net,tt)            = sum(time_map(tt,ttX), link_capX(net,ttX) * overlap(tt,ttX));
link_efficiency2(net,tt)     =   sum(time_map(tt,ttX), link_efficiencyX(net,ttX) * overlap(tt,ttX))
                               / sum(time_map(tt,ttX), overlap(tt,ttX));
demand2(rr,tt)               = sum(time_map(tt,ttX), demandX(rr,ttX) * overlap(tt,ttX));
avail2(rr,p,tt)              = sum(time_map(tt,ttX), availX(ttX,rr,p)) / sum(time_map(tt,ttX), 1);
$ifi not %METHOD%==spExplicitDE type_mult(type)= 1;

* projection on parameters with rearranged index sets
option plant_cap<plant_cap2, link_cap<link_cap2, link_efficiency<link_efficiency2, avail<avail2;
option demand<demand2;

* *** ERRORCHECKS
set err01(rr,p) 'consistent region plant type data';
err01(rr,p) = rp(rr,p) xor sum(ptype(rr,p,type), 1);
abort$card(err01) err01, rp, ptype;


* ##############################################################################
* if method=pips prepare regional blocks (rb) and time blocks (tb) and map them to the main blocks (bb)
$IFTHENi.method %METHOD%==pips
   Scalar BLOCK0 /1/, BLOCKN;
   BLOCKN = BLOCK0 + card(rr);
   file fopt / "%gams.optdir%convertd.opt" /;
$  if not set rbsize $set  rbsize  %NBREGIONS%
$  if not set tbsize $set  tbsize  24
$  eval rblock  ceil(card(rr)/%rbsize%)
$  eval tblock  ceil(card(tt)/%tbsize%)
$  log %rblock%
$  log %tblock%
   set rb regional block / rb1*rb%rblock% /
       tb time block     / tb1*tb%tblock% /;
$  eval bmax %rblock%*%tblock%
$  log %bmax%
   set
       bb                / b1*b%bmax% /
       b(bb)
       rcnt              / rcnt1*rcnt%rbsize% /
       tcnt              / tcnt1*tcnt%tbsize% /
       bmap(bb,rb,tb)    / #bb:(#rb.#tb) /
       rmap(rr,rb,rcnt)  / #rr:(#rb.#rcnt) /
*       rmap(rr,rb,rcnt)  / r1.rb1.rcnt1, r4.rb1.rcnt2, r6.rb1.rcnt3
*                           r2.rb2.rcnt1, r3.rb2.rcnt2, r5.rb2.rcnt3
*                           r7.rb3.rcnt1   /
       tmap(tt,tb,tcnt)  / #tt:(#tb.#tcnt) /
       xmap(bb,rb,tb,rb,rr,tb,tt)
       brtmap(bb,rr,tt)
       btrmap(bb,tt,rr)
       rtbmap(rr,tt,bb)
       trbmap(tt,rr,bb)
       rbrmap(rb,rr)
       tbtmap(tb,tt)
       ttbmap(tt,tb)
       tsubseq(tt1,tt2)
   ;
   alias(rtbmap,rtbmap2), (trbmap,trbmap2);
   tsubseq(tt,tt++1) = yes;
   option rbrmap<rmap, tbtmap<tmap;
   xmap(bmap(bb,rb,tb),rbrmap(rb,rr),tbtmap(tb,tt)) = yes;
   option brtmap<xmap, btrmap<brtmap, trbmap<btrmap, rtbmap<brtmap, ttbmap<tbtmap;

   parameter rt_stage(rr,tt)
             tr_stage(tt,rr)
             rrt_stage(rr1,rr2,tt)
             trr_stage(tt,rr1,rr2)
             rtt_stage(rr,tt1,tt2)
             rt2_stage(rr1,rr2,tt1,tt2)
   ;

   rt_stage(rr,tt)                           = sum(rtbmap(rr,tt,bb),ord(bb)) + 1;
   rrt_stage(net(rr1,rr2),tt)                = sum((rtbmap(rr1,tt,bb),rtbmap2(rr2,tt,bb)),  ord(bb)) + 1;
   rtt_stage(rr,tsubseq(tt1,tt2))            = sum((rtbmap(rr,tt1,bb),rtbmap2(rr,tt2,bb)),  ord(bb)) + 1;
   option tr_stage<rt_stage, trr_stage<rrt_stage;
$else.method
  set bb    /b0/
      b(bb) /b0/
      btrmap(bb,tt,rr)
  ;
  btrmap(bb,tt,rr) = yes;
$endif.method
alias (bb,bb1,bb2);

* ##############################################################################

Scalar genBlock indicator for generating a block in case of method=PIPS /0/;


Positive variables
    POWER(tt,rr,p)           'power production                   [GWh]'
    FLOW(tt,rr1,rr2)         'power flow                         [GWh]'
    LINK_ADD_CAP(rr1,rr2)    'arc capacity expansion             [GWh]'
    SLACK(tt,rr)             'uncovered demand                   [GWh]'
    STORAGE_LEVEL(tt,rr,s)   'storage level                      [GWh]'
    STORAGE_INFLOW(tt,rr,s)  'power entering storage             [GWh]'
    STORAGE_OUTFLOW(tt,rr,s) 'power taken from storage           [GWh]'
    PLANT_ADD_CAP(rr,p)      'plant capacity expansion            [GW]'
    STORAGE_ADD_CAP(rr,s)    'storage capacity expansion         [GWh]'
Variable
    DUMMY
    ROBJ(rr)                 'total region cost                 [MEUR]'
    BOBJ(bb)                 'total block cost                  [MEUR]'
    OBJ                      'total                             [MEUR]'
;

* bounds
LINK_ADD_CAP.up(net)      = link_max_add_cap(net)  ;
PLANT_ADD_CAP.up(rp)      = plant_max_add_cap(rp)  ;
STORAGE_ADD_CAP.up(rs)    = storage_max_add_cap(rs);
STORAGE_INFLOW.up(tt,rs)  = storage_max_in(rs)     ;
STORAGE_OUTFLOW.up(tt,rs) = storage_max_out(rs)    ;

equations
    eq_bobj(bb)                     'total cost in block                  '
    eq_power_balance(tt,rr)         'power balance                        '
    eq_plant_capacity(tt,rr,p)      'respect plant capacity               '
    eq_total_plant_capacity(rr,p)   'respect total plant capacity         '
    eq_storage_balance(tt,rr,s)     'storage balance                      '
    eq_storage_capacity(tt,rr,s)    'respect storage capacity             '
    eq_emission_cap(e)              'respect emission cap                 '
    eq_link_capacity(tt,rr1,rr2)    'respect link capacity                '
    eq_obj                          'total cost                           '
;

* subsets driving the equations
set s_eq_bobj(bb)
    s_eq_power_balance(tt,rr)
    s_eq_plant_capacity(tt,rr,p)
    s_eq_total_plant_capacity(rr,p)
    s_eq_storage_balance(tt,rr,s)
    s_eq_storage_capacity(tt,rr,s)
    s_eq_emission_cap(e)
    s_eq_link_capacity(tt,rr1,rr2)
;

eq_obj..
    OBJ =e= sum(b, BOBJ(b))
          + sum(net(r1,r2), LINK_ADD_CAP(net)   * cost_link_add(net))
          + sum(rp(r,p),    PLANT_ADD_CAP(rp)   * cost_plant_add(rp))
          + sum(rs(r,s),    STORAGE_ADD_CAP(rs) * cost_storage_add(rs))
$iftheni.method %METHOD%==pips
          // list of linking variables that actually do not appear in the objective
          + eps * sum((tt,net)$[FLOW.stage(tt,net)=1], FLOW(tt,net))
          + eps * sum(net$[LINK_ADD_CAP.stage(net)=1], LINK_ADD_CAP(net))
          + eps * sum(rp(rr,p)$[PLANT_ADD_CAP.stage(rp)=1], PLANT_ADD_CAP(rp))
          + eps * sum(rs(rr,s)$[STORAGE_ADD_CAP.stage(rs)=1], STORAGE_ADD_CAP(rs))
$endif.method
;

eq_bobj(s_eq_bobj(b))..
$ifthene.scaling NOT %SCALING% == 1
    BOBJ(b) =e= sum(btrmap(b,t,r), sum((ptype(rp(r,p),type)), POWER(t,rp) * cost_power_generation(rp) * type_mult(type))
                                 + SLACK(t,r) * cost_unserved_demand(t)
                                 + sum((rp(r,p),e), POWER(t,rp) * plant_emission(rp,e) * cost_emission(e)) );
$else.scaling
   [BOBJ(b) - sum(btrmap(b,t,r), sum((ptype(rp(r,p),type)), POWER(t,rp) * cost_power_generation(rp) * type_mult(type))
                               + SLACK(t,r) * cost_unserved_demand(t)
                               + sum((rp(r,p),e), POWER(t,rp) * plant_emission(rp,e) * cost_emission(e)) )]
   / [10*smin((btrmap(b,t,r),ptype(rp(r,p),type),e)$cost_power_generation(rp),  cost_power_generation(rp)*type_mult(type) + plant_emission(rp,e) * cost_emission(e))]
   =e= 0;
$endif.scaling

eq_power_balance(s_eq_power_balance(t,r))..
        sum(rp(r,p),    POWER(t,rp))
      + sum(net(rr2,r), FLOW(t,net) * link_efficiency(t,net))
      - sum(net(r,rr2), FLOW(t,net))
      + sum(rs(r,s),    STORAGE_OUTFLOW(t,rs) - STORAGE_INFLOW(t,rs))
      + SLACK(t,r)
    =g= demand(t,r);

eq_plant_capacity(s_eq_plant_capacity(t,rp(r,p)))..
$ifthene.scaling NOT %SCALING% == 1
    POWER(t,rp) =l= (plant_cap(t,rp) + PLANT_ADD_CAP(rp)*%RESOLUTION%) * avail(t,rp) ;
$else.scaling
    [POWER(t,rp) - (plant_cap(t,rp) + PLANT_ADD_CAP(rp)*%RESOLUTION%) * avail(t,rp)] / [1$(avail(t,rp)<1e-6) + ((%RESOLUTION%)*avail(t,rp)*10)$(avail(t,rp)>=1e-6)] =l= 0;
$endif.scaling

eq_total_plant_capacity(s_eq_total_plant_capacity(rp(r,p)))..
    sum(t, POWER(t,rp)) =l= total_plant_cap(rp);

eq_storage_balance(s_eq_storage_balance(t(tt),rs(r,s)))..
    STORAGE_LEVEL(t,rs) =e= STORAGE_LEVEL(tt--1,rs)$t2(tt--1) * storage_efficiency(rs)
                          + STORAGE_INFLOW(t,rs)    * storage_efficiency_in(rs)
                          - STORAGE_OUTFLOW(t,rs)   / storage_efficiency_out(rs);

eq_storage_capacity(s_eq_storage_capacity(t,rs(r,s)))..
    STORAGE_LEVEL(t,rs) =l= storage_cap(rs) + STORAGE_ADD_CAP(rs);

eq_emission_cap(s_eq_emission_cap(e))..
$ifthene.scaling NOT %SCALING% == 1
    sum((rp(r,p),t), POWER(t,rp) * plant_emission(rp,e)) =l= total_emission_cap(e);
$else.scaling
    [sum((rp(r,p),t), POWER(t,rp) * plant_emission(rp,e)) - total_emission_cap(e)]
    / [10*smin((rp(r,p),t)$plant_emission(rp,e), plant_emission(rp,e))]
    =l= 0;
$endif.scaling





eq_link_capacity(s_eq_link_capacity(t,net))..
    FLOW(t,net) =l= link_cap(t,net) + LINK_ADD_CAP(net) * %RESOLUTION%;

model simple / all /;

$ifthene.noslack %NOSLACK%==1
  SLACK.fx(tt,rr)  = 0;
  simple.holdfixed = 1;
$endif.noslack

$iftheni.method %METHOD%==pips
   option lp=convertd;
   simple.optfile=1;
   simple.solvelink=2;
   simple.limrow=0;
   simple.limcol=0;
   simple.solprint=0;

   b(bb) = yes;
   t(tt) = yes;
   r(rr) = yes;
   //Variables
   POWER.stage(t,rp(r,p))            = sum(trbmap(t,r,b(bb)), ord(bb)) + 1;
   FLOW.stage(t,net(r1,r2))          = sum((trbmap(t,r1,b(bb)),trbmap2(t,r2,b)), ord(bb)) + 1;
   LINK_ADD_CAP.stage(net(r1,r2))    = smin(rtbmap(r1,tt,b(bb)), ord(bb))$[card(tb)<=1 and sum((rtbmap(r1,tt,b),rtbmap2(r2,tt,b)), 1)] + 1;
   SLACK.stage(t,r)                  = sum(trbmap(t,r,bb), ord(bb)) + 1;
   STORAGE_LEVEL.stage(t,rs(r,s))    = sum(trbmap(t,r,bb), ord(bb)) + 1;
   STORAGE_INFLOW.stage(t,rs(r,s))   = sum(trbmap(t,r,bb), ord(bb)) + 1;
   STORAGE_OUTFLOW.stage(t,rs(r,s))  = sum(trbmap(t,r,bb), ord(bb)) + 1;
   PLANT_ADD_CAP.stage(rp(r,p))      = smin(rtbmap(r,tt,b(bb)), ord(bb))$[card(tb)<=1] + 1; // if no time decomposition put in region block
   STORAGE_ADD_CAP.stage(rs(r,s))    = smin(rtbmap(r,tt,b(bb)), ord(bb))$[card(tb)<=1] + 1; // if no time decomposition put in region block
   BOBJ.stage(b(bb))                 = ord(bb) + 1;
   //Equations
   eq_bobj.stage(b(bb))                   = ord(bb) + 1;
   eq_power_balance.stage(t,r)            = sum(trbmap(t,r,b(bb)), ord(bb)) + 1;
   eq_plant_capacity.stage(t,rp(r,p))     = sum(trbmap(t,r,b(bb)), ord(bb)) + 1;
   eq_total_plant_capacity.stage(rp(r,p)) =   smin(rtbmap(r,tt,b(bb)), ord(bb) + 1)$[card(tb)<=1] // no time decomposition --> region block
                                            + (card(bb) + 2)$[card(tb)>=2];                       // time decomposition    --> coupling constraint
   eq_storage_balance.stage(t,rs(r,s))$[    sum((ttbmap(t,tb),tsubseq(tt,t))$tbtmap(tb,tt),1)] = sum(trbmap(t,r,b(bb)), ord(bb)) + 1;
   eq_storage_balance.stage(t,rs(r,s))$[not sum((ttbmap(t,tb),tsubseq(tt,t))$tbtmap(tb,tt),1)] = card(bb) + 2;
   eq_storage_capacity.stage(t,rs(r,s))   =  sum(trbmap(t,r,b(bb)), ord(bb)) + 1;
   eq_emission_cap.stage(e)               = card(bb) + 2;
   eq_link_capacity.stage(t,net(r1,r2))   = sum((trbmap(t,r1,b(bb)),trbmap2(t,r2,b)), ord(bb)) + 1;
   eq_obj.stage = card(bb) + 2;
   //equation controlling subsets
   s_eq_bobj(b)                       = yes;
   s_eq_power_balance(t,r)            = yes;
   s_eq_plant_capacity(t,rp(r,p))     = yes;
   s_eq_total_plant_capacity(rp(r,p)) = yes;
   s_eq_storage_balance(t,rs(r,s))    = yes;
   s_eq_storage_capacity(t,rs(r,s))   = yes;
   s_eq_emission_cap(e)               = yes;
   s_eq_link_capacity(t,net(r1,r2))   = yes;

   putclose fopt 'jacobian allblocksPips.gdx';
*  / 'dictmap dmallblocksPips.gdx';
   solve simple min OBJ use lp;

$ELSE.method
*  Standard LP
   t(tt) = yes;
   r(rr) = yes;

   s_eq_bobj(b)                       = yes;
   s_eq_power_balance(t,r)            = yes;
   s_eq_plant_capacity(t,rp(r,p))     = yes;
   s_eq_total_plant_capacity(rp(r,p)) = yes;
   s_eq_storage_balance(t,rs(r,s))    = yes;
   s_eq_storage_capacity(t,rs(r,s))   = yes;
   s_eq_emission_cap(e)               = yes;
   s_eq_link_capacity(t,net(r1,r2))   = yes;

   simple.optfile = 1;
$onecho > cplex.opt
lpmethod 4
solutiontype 2
preind 0
$offecho

   solve simple min OBJ use lp;

$ENDIF.method


