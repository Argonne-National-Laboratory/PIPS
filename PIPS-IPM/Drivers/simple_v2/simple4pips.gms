$ eolcom //

$if  not set FROM                 $set FROM              0
$if  not set TO                   $set TO                1
$if  not set RESOLUTION           $set RESOLUTION        1
$if  not set NBREGIONS            $set NBREGIONS         4
$if  not set NBSHIFTS             $set NBSHIFTS          0
$if  not set SHIFTSTEPSIZE        $set SHIFTSTEPSIZE     2
$if  not set METHOD               $set METHOD            standard_lp
$if  not set BLOCK                $set BLOCK            -2
$if  not set SCALING              $set SCALING           0
$if  not set NOSLACK              $set NOSLACK           0
$if  not set KEEPVENAMES          $set KEEPVENAMES       0
$if  not set KEEPUELS             $set KEEPUELS          1
$if  not set SUPPRESSDM           $set SUPPRESSDM        1
$ifi not set SLICE                $set SLICE             0
$ifi not set RUNPIPS              $set RUNPIPS           0

$ifi not %METHOD%==PIPS           $set SLICE             0
$ife %BLOCK%<=-1                  $set SLICE             0

$setglobal GDXSUFFIX
$ife %KEEPVENAMES%==0 $setglobal GDXSUFFIX %GDXSUFFIX%_noVEnames
$ife %KEEPUELS%==0    $setglobal GDXSUFFIX %GDXSUFFIX%_noUELs

$ife %FROM%>%TO%         $abort 'FROM > TO'
$ife %RESOLUTION%<0      $abort 'Negative RESOLUTION forbidden'

$eval NBTIMESTEPS ceil((%TO%-%FROM%)*365*24/%RESOLUTION%)
$setnames "%gams.input%" SIMPLEDIR fn fe

*   basic sets
Set rr                          'regions                                       '
    p                           'plants                                        '
    s                           'storages                                      '
    tt                          'time steps                                    ' / t1*t%NBTIMESTEPS% /
    ttX                         'time steps in input data                      '
    e                           'emissions                                     '
    type                        'plant type                                    '
    bb                          'blocks                                        '
;
singleton set ttmin(tt);
ttmin(tt) = ord(tt) = 1;
Alias(rr,rr1,rr2), (tt,tt1,tt2);

*   Subsets and mappings
set r(rr)                       'active region                                 '
    t(tt)                       'subset of active time steps                   '
    t_strg_linking(tt)          'time steps for which strg balance is linking  '
    rp(rr,p)                    'region to plant mapping                       '
    ptype(rr,p,type)            'plant type mapping                            '
    rs(rr,s)                    'region to storage mapping                     '
    net(rr1,rr2)                'transmission links                            '
    time_map(tt,ttX)            'mapping of model to input time steps          '
    tX(ttX)                     'dynamic subset of data time steps             '
    tlast(tt)                   'last time step                                '
    tXlast(ttX)                 'last hour in base data                        '
$ifthene.shft %NBSHIFTS%>0
    shiftclass                  'shift classes                                 ' / shift1*shift%NBSHIFTS% /
$else.shft
$   onempty
    shiftclass(*)               'shift classes                                 ' / /
$   offempty
$endif.shft
    t_shift_linking(tt,shiftclass) 'timesteps where ecar balance constraint is linking'
;
Alias(r,r1,r2),(t,t1,t2);

Parameter
    plant_capX(rr,p,ttX)         'plant capacity                           [GW]'
    cost_unserved_demandX(ttX)   'price for unserved demand          [MEUR/GWh]'
    link_capX(rr1,rr2,ttX)       'transmission link capacity per hour     [GWh]'
    link_efficiencyX(rr1,rr2,ttX)'transmission link efficiency factor          '
    demandX(rr,ttX)              'demand for region per hour              [GWh]'
    demandEcarX(rr,ttX)          'ecar demand for region per hour         [GWh]'
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
    demandEcar(tt,rr)            'ecar demand for region per time step    [GWh]'
    demandEcar2(rr,tt)           'ecar demand for region per time step    [GWh]'
    plant_max_add_cap(rr,p)      'max additional plant capacity            [GW]'
    storage_max_add_cap(rr,s)    'max additional storage capacity         [GWh]'
    link_max_add_cap(rr1,rr2)    'max additional arc capacity             [GWh]'
    cost_plant_add(rr,p)         'cost for additional plant capacity  [MEUR/GW]'
    cost_storage_add(rr,s)       'cost for additional storage cap    [MEUR/GWh]'
    cost_link_add(rr1,rr2)       'cost for additional link capacity  [MEUR/GWh]'
    type_mult(type)              'generation cost multiplier                   '
    availX(ttX,rr,p)             'plant availibility vector                    '
    avail(tt,rr,p)               'plant availibility vector                    '
    avail2(rr,p,tt)              'plant availibility vector                    '
    shiftSteps(shiftclass)       'number of time steps per shiftclass          '
;
shiftSteps(shiftclass) = ord(shiftclass)*%SHIFTSTEPSIZE%;

scalar
    tmpStart                     'helper storing time step start times         '
    tmpEnd                       'helper storing time step end times           '
    shiftlimit                   'maximum share of ecar demand that can be shifted' / 0.5 /           
    shiftCost                    'cost for shifting demand satisfaction to later time step' / 0.01 /
;

$iftheni.method %METHOD%==PIPS
*  compute time, region and total number of blocks
$  if not set rbsize $set  rbsize  %NBREGIONS%
$  if not set tbsize $set  tbsize  24
$  eval rblock  ceil(%NBREGIONS%/%rbsize%)
$  eval tblock  ceil(%NBTIMESTEPS%/%tbsize%)
$  eval bmax %rblock%*%tblock%
$  log #region blocks: %rblock%
$  log #time blocks:   %tblock%
$  log #blocks:        %bmax%
$endif.method


$ifthene.readGdx %SLICE%==2
*  if data slicing level 2 is enabled, get the data from gdx (must be created in advance)
$  eval txFirst floor((%RESOLUTION%*%TBSIZE%)*(%BLOCK%-1)+1)
$  eval txLast   ceil((%RESOLUTION%*%TBSIZE%)*(%BLOCK%))
$  log #txFirst:       %txFirst%
$  log #txLast:        %txLast%
$  ifthene.blk %BLOCK%>=1
     set ttX / tX%tXfirst%*tX%tXlast% /;
$  else.blk
$    onempty
     set ttX / /;
$    offempty
$  endif.blk
$  if not exist simple_data_%NBREGIONS%.gdx $abort SLICE=2 requires previous generation of simple_data_%NBREGIONS%.gdx
$  gdxin  "simple_data_%NBREGIONS%.gdx"
$  load startX endX
$  load rr
$  load p s e rp type ptype rs net
$  load plant_capX cost_unserved_demandX link_capX link_efficiencyX demandX demandEcarX
$  load yearly_plant_cap plant_emission cost_power_generation yearly_emission_cap
$  load cost_emission storage_cap storage_efficiency storage_efficiency_in
$  load storage_efficiency_out storage_max_in storage_max_out
$  load plant_max_add_cap storage_max_add_cap link_max_add_cap cost_plant_add
$  load cost_storage_add cost_link_add availX
$  gdxin
$else.readGdx
*  else run data generator as part of simple4pips
$  include "%SIMPLEDIR%simple_data_gen.gms"
$endif.readGdx

* ##############################################################################
$IFTHENi.method %METHOD%==PIPS
*  create block sets and compute mappings
   Scalar BLOCK0 / 1 /, BLOCKN;
   BLOCKN = BLOCK0 + card(rr);
   file fopt / "%gams.optdir%convertd.opt" /;
   set rb regional block / rb1*rb%rblock% /
       tb time block     / tb1*tb%tblock% /
       bb                / b1*b%bmax% /
       b(bb)
       rcnt              / rcnt1*rcnt%rbsize% /
       tcnt              / tcnt1*tcnt%tbsize% /
       bmap(bb,rb,tb)    / #bb:(#rb.#tb) /
       rmap(rr,rb,rcnt)  / #rr:(#rb.#rcnt) /
       tmap(tt,tb,tcnt)  / #tt:(#tb.#tcnt) /
       bf               'block files       ' / block0*block%bmax% /
       actbf(bf)        'active block files'
       xmap(bb,rb,tb,rb,rr,tb,tt)
       brtmap(bb,rr,tt)
       btrmap(bb,tt,rr)
       rtbmap(rr,tt,bb)
       trbmap(tt,rr,bb)
       rbrmap(rb,rr)
       rrbmap(rr,rb)
       tbtmap(tb,tt)
       ttbmap(tt,tb)
       tsubseq(tt1,tt2)
       samerb(rr,rr)
   ;
   alias(rtbmap,rtbmap2);
   tsubseq(tt,tt++1) = yes;
   option rbrmap<rmap, rrbmap<rmap, tbtmap<tmap;
$  ifthene.slice %SLICE%>=1
*    if data slicing level 1 enabled, generate all data but then use dynamic subset of t with relevent time steps
$    ifthene.blk %BLOCK%>=1
       t(tt) = sum(tbtmap('tb%BLOCK%',tt), 1);
$    else.blk
       t(tt) = no;
$    endif.blk
$  else.slice
*    else, use all time steps
     t(tt) = yes;
$  endif.slice
   xmap(bmap(bb,rb,tb),rbrmap(rb,rr),tbtmap(tb,t)) = yes;
   option brtmap<xmap, btrmap<brtmap, trbmap<btrmap, rtbmap<brtmap, ttbmap<tbtmap;
   option clear=xmap;
   t_strg_linking(tt1)  = not sum((ttbmap(tt1,tb),tsubseq(tt2,tt1))$tbtmap(tb,tt2),1); //identify time steps where storage balance is linking constraint
   t_shift_linking(tt1,shiftclass) = not sum(ttbmap(tt1,tb)$tbtmap(tb,tt1--shiftSteps(shiftclass)), 1);
   samerb(net(rr1,rr2)) = sum(rb$(rrbmap(rr1,rb) and rrbmap(rr2,rb)), 1);
   //compute stage annotations (once) for certain important combinations of index sets
   parameter rt_stage(rr,tt)
             tr_stage(tt,rr)
             rrt_stage(rr1,rr2,tt)
             trr_stage(tt,rr1,rr2)
   ;
   rt_stage(rr,t(tt))                              = sum(rtbmap(rr,tt,bb),ord(bb)) + 1;
   rrt_stage(net(rr1,rr2),t(tt))$(not samerb(net)) = 1;                //regions in different rb
   rrt_stage(net(samerb(rr1,rr2)),t(tt))           = rt_stage(rr1,tt); //regions in same rb
   option tr_stage<rt_stage, trr_stage<rrt_stage;
$else.method
*  if PIPS not activated use "one block" containing all regions and tinme steps
   set bb    /b0/
       b(bb) /b0/
       btrmap(bb,tt,rr)
   ;
   btrmap(bb,tt,rr) = yes;
   t(tt) = yes;
$endif.method

* ##############################################################################


tlast(tt)  = tt.last;
* compute mapping and overlap of model time steps to input time steps
start(tt)   = %FROM%*365*24 + (ord(tt)-1) * %RESOLUTION%;
end(tt)     = %FROM%*365*24 +  ord(tt)    * %RESOLUTION%;
end(tlast)  = %TO%*365*24;
tX(ttX)     = no;

loop(ttX$(not card(tX)), tX(ttX) = endX(ttX) > smin(t(tt), start(tt)););
tmpStart = sum(tX, startX(tX));
tmpEnd   = sum(tX, endX(tX));

loop(t(tt),
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
$ifthene not %SLICE%>=1
abort$(not sum(time_map(tlast,tXlast),1)) 'Last time steps of base and model data do not map.';
$endif
overlap(time_map(tt,ttX)) = min(1, abs(min(0, max(end(tt),endX(ttX)) - min(start(tt),startX(ttX)) - (%RESOLUTION%+1))));
overlap(time_map(tlast,ttX)) = min(end(tlast),endX(ttX)) - max(start(tlast),startX(ttX));

* Compute parameters according to time span and resolution
plant_cap2(rp(rr,p),t(tt))   = sum(time_map(tt,ttX), plant_capX(rp,ttX) * overlap(tt,ttX));
total_plant_cap(rp(rr,p))    = (%TO% - %FROM%) * yearly_plant_cap(rp);
total_emission_cap(e)        = (%TO% - %FROM%) * yearly_emission_cap(e);
storage_efficiency(rs(rr,s)) = rPower(storage_efficiency(rs),%RESOLUTION%);
storage_max_in(rs(rr,s))     = %RESOLUTION% * storage_max_in(rs);
storage_max_out(rs(rr,s))    = %RESOLUTION% * storage_max_out(rs);
cost_unserved_demand(t(tt))  = smax(time_map(tt,ttX), cost_unserved_demandX(ttX));
link_cap2(net,t(tt))         = sum(time_map(tt,ttX), link_capX(net,ttX) * overlap(tt,ttX));
link_efficiency2(net,t(tt))  =   sum(time_map(tt,ttX), link_efficiencyX(net,ttX) * overlap(tt,ttX))
                               / sum(time_map(tt,ttX), overlap(tt,ttX));
demand2(rr,t(tt))            = sum(time_map(tt,ttX), demandX(rr,ttX) * overlap(tt,ttX));
demandEcar2(rr,t(tt))        = sum(time_map(tt,ttX), demandEcarX(rr,ttX) * overlap(tt,ttX));
avail2(rp(rr,p),t(tt))       = sum(time_map(tt,ttX), availX(ttX,rp)) / sum(time_map(tt,ttX), 1);
$ifi not %METHOD%==spExplicitDE type_mult(type)= 1;

option clear=availX, clear=plant_capX, clear=link_capX, clear=link_efficiencyX, clear=demandX, clear=demandEcarX; //clear hourly input data
* projection on parameters with rearranged index sets
option plant_cap<plant_cap2, link_cap<link_cap2, link_efficiency<link_efficiency2, avail<avail2;
option demand<demand2, demandEcar<demandEcar2;
option clear=link_cap2, clear=link_efficiency2, clear=demand2, clear=demandEcar2, clear=avail2; //clear helpers used for speedy computations

* ERRORCHECKS
set err01(rr,p) 'inconsistent region plant type data';
err01(rr,p) = rp(rr,p) xor sum(ptype(rr,p,type), 1);
abort$card(err01) err01, rp, ptype;

Scalar
    equStage1 'dynamic indicator for linking equation stage' / -1 /
;

*macro to add DUMMY Var to an equation if it is linking
$macro linkEqEpsDummy(equName,arglist)  (eps*DUMMY)$(equName.stage&&arglist =equStage1)

*macros to generate equation/variable only if index sets are in dynamic domain
$macro genEqu(equName,domain,arglist)   equName&&domain$(gen_&&equName&&arglist)
$macro genVar(varName,arglist)          varName&&arglist$(gen_&&varName&&arglist)

*macro with list of (potentially) linking variables
alias(tt,ttM), (net,netM), (rr,rrM), (p,pM), (rp,rpM), (s,sM), (rs,rsM)
$macro epsLinkVar   eps * sum((ttM,netM)$[lvar_FLOW(ttM,netM)], FLOW(ttM,netM))              \
                  + eps * sum(netM$[lvar_LINK_ADD_CAP(netM)], LINK_ADD_CAP(netM))            \
                  + eps * sum(rpM(rrM,pM)$[lvar_PLANT_ADD_CAP(rpM)], PLANT_ADD_CAP(rpM))     \
                  + eps * sum(rsM(rrM,sM)$[lvar_STORAGE_ADD_CAP(rsM)], STORAGE_ADD_CAP(rsM)) \


Positive variables
    POWER(tt,rr,p)                 'power production                               [GWh]'
    ECAR_BALANCE(tt,rr,shiftclass) 'additional ecar demand from earlier time steps [GWh]'
    ECAR_DELAY(tt,rr,shiftclass)   'ecar demand delayed to later time steps        [GWh]'
    FLOW(tt,rr1,rr2)               'power flow                                     [GWh]'
    LINK_ADD_CAP(rr1,rr2)          'arc capacity expansion                         [GWh]'
    SLACK(tt,rr)                   'uncovered demand                               [GWh]'
    STORAGE_LEVEL(tt,rr,s)         'storage level                                  [GWh]'
    STORAGE_INFLOW(tt,rr,s)        'power entering storage                         [GWh]'
    STORAGE_OUTFLOW(tt,rr,s)       'power taken from storage                       [GWh]'
    PLANT_ADD_CAP(rr,p)            'plant capacity expansion                        [GW]'
    STORAGE_ADD_CAP(rr,s)          'storage capacity expansion                     [GWh]'
Variable                                                                        
    DUMMY                          'dummy linking variable                              '
    BOBJ(bb)                       'total block cost                              [MEUR]'
    OBJ                            'total                                         [MEUR]'
;

* bounds
LINK_ADD_CAP.up(net)         = link_max_add_cap(net)  ;
PLANT_ADD_CAP.up(rp)         = plant_max_add_cap(rp)  ;
STORAGE_ADD_CAP.up(rs)       = storage_max_add_cap(rs);
STORAGE_INFLOW.up(t(tt),rs)  = storage_max_in(rs)     ;
STORAGE_OUTFLOW.up(t(tt),rs) = storage_max_out(rs)    ;

equations
    eq_bobj(bb)                       'total cost in block              '
    eq_power_balance(tt,rr)           'power balance                    '
    eq_ecar_balance(tt,rr,shiftclass) 'compute additional ecar demand due to shifting from earlier time steps'
    eq_ecar_shift_limit(tt,rr)        'limit shiftting ecar demand to later time steps'
    eq_plant_capacity(tt,rr,p)        'respect plant capacity           '
    eq_total_plant_capacity(rr,p)     'respect total plant capacity     '
    eq_storage_balance(tt,rr,s)       'storage balance                  '
    eq_storage_capacity(tt,rr,s)      'respect storage capacity         '
    eq_emission_cap(e)                'respect emission cap             '
    eq_link_capacity(tt,rr1,rr2)      'respect link capacity            '
    eq_obj                            'total cost                       '
;

*dynamic domain sets for variables, equations and linking variables
set gen_POWER(tt,rr,p)                       'dynamic variable domain set'
    gen_ECAR_BALANCE(tt,rr,shiftclass)       'dynamic variable domain set'
    gen_ECAR_DELAY(tt,rr,shiftclass)         'dynamic variable domain set'
    gen_FLOW(tt,rr1,rr2)                     'dynamic variable domain set'
    gen_LINK_ADD_CAP(rr1,rr2)                'dynamic variable domain set'
    gen_SLACK(tt,rr)                         'dynamic variable domain set'
    gen_STORAGE_LEVEL(tt,rr,s)               'dynamic variable domain set'
    gen_STORAGE_INFLOW(tt,rr,s)              'dynamic variable domain set'
    gen_STORAGE_OUTFLOW(tt,rr,s)             'dynamic variable domain set'
    gen_PLANT_ADD_CAP(rr,p)                  'dynamic variable domain set'
    gen_STORAGE_ADD_CAP(rr,s)                'dynamic variable domain set'
    gen_ROBJ(rr)                             'dynamic variable domain set'
    gen_BOBJ(bb)                             'dynamic variable domain set'
    gen_eq_bobj(bb)                          'dynamic equation domain set'
    gen_eq_power_balance(tt,rr)              'dynamic equation domain set'
    gen_eq_ecar_balance(tt,rr,shiftclass)    'dynamic equation domain set'
    gen_eq_ecar_shift_limit(tt,rr)           'dynamic equation domain set'
    gen_eq_plant_capacity(tt,rr,p)           'dynamic equation domain set'
    gen_eq_total_plant_capacity(rr,p)        'dynamic equation domain set'
    gen_eq_storage_balance(tt,rr,s)          'dynamic equation domain set'
    gen_eq_storage_capacity(tt,rr,s)         'dynamic equation domain set'
    gen_eq_emission_cap(e)                   'dynamic equation domain set'
    gen_eq_link_capacity(tt,rr1,rr2)         'dynamic equation domain set'
    lvar_FLOW(tt,rr,rr)                      'dynamic linking variable domain set'   
    lvar_LINK_ADD_CAP(rr,rrM)                'dynamic linking variable domain set'
    lvar_PLANT_ADD_CAP(rr,p)                 'dynamic linking variable domain set'
    lvar_STORAGE_ADD_CAP(rr,s)               'dynamic linking variable domain set'
;

eq_obj..
    OBJ + linkEqEpsDummy(eq_obj,'') =e= sum(b, genVar(BOBJ,(b)))
                                      + sum(net(r1,r2), genVar(LINK_ADD_CAP,(net))   * cost_link_add(net))
                                      + sum(rp(r,p),    genVar(PLANT_ADD_CAP,(rp))   * cost_plant_add(rp))
                                      + sum(rs(r,s),    genVar(STORAGE_ADD_CAP,(rs)) * cost_storage_add(rs))
                                      + epsLinkVar ;

genEqu(eq_bobj,(b),(b))..
$ifthene.scaling NOT %SCALING% == 1
    genVar(BOBJ,(b)) =e= sum(btrmap(b,t,r), sum((ptype(rp(r,p),type)), genVar(POWER,(t,rp)) * cost_power_generation(rp) * type_mult(type))
                                               + genVar(SLACK,(t,r)) * cost_unserved_demand(t)
                                               + sum((rp(r,p),e), genVar(POWER,(t,rp)) * plant_emission(rp,e) * cost_emission(e))
                                               + shiftCost*sum(shiftclass, ECAR_DELAY(t,r,shiftclass)) );
$else.scaling
     [genVar(BOBJ,(b)) - sum(btrmap(b,t,r), sum((ptype(rp(r,p),type)), genVar(POWER,(t,rp)) * cost_power_generation(rp) * type_mult(type))
                       + genVar(SLACK,(t,r)) * cost_unserved_demand(t)
                       + sum((rp(r,p),e), genVar(POWER,(t,rp)) * plant_emission(rp,e) * cost_emission(e)) )]
   / [10*smin((btrmap(b,t,r),ptype(rp(r,p),type),e)$cost_power_generation(rp),  cost_power_generation(rp)*type_mult(type) + plant_emission(rp,e) * cost_emission(e))]
   =e= 0;
$endif.scaling

genEqu(eq_power_balance,(t,r),(t,r))..
        sum(rp(r,p),    genVar(POWER,(t,rp)))
      + sum(net(rr2,r), genVar(FLOW,(t,net)) * link_efficiency(t,net))
      - sum(net(r,rr2), genVar(FLOW,(t,net)))
      + sum(rs(r,s),    genVar(STORAGE_OUTFLOW,(t,rs)) - genVar(STORAGE_INFLOW,(t,rs)))
      + genVar(SLACK,(t,r))
      + linkEqEpsDummy(eq_power_balance,(t,r))
    =g= demand(t,r)
      + demandEcar(t,r) + sum(shiftclass, ECAR_BALANCE(t,r,shiftclass) - ECAR_DELAY(t,r,shiftclass)) ;


genEqu(eq_ecar_balance,(tt,r,shiftclass),(tt,r,shiftclass))..
    genVar(ECAR_BALANCE,(tt,r,shiftclass)) + linkEqEpsDummy(eq_ecar_balance,(tt,r,shiftclass)) =e= genVar(ECAR_DELAY,(tt--shiftSteps(shiftclass),r,shiftclass));

genEqu(eq_ecar_shift_limit,(t,r),(t,r))..
    sum(shiftclass, genVar(ECAR_DELAY,(t,r,shiftclass))) =l= shiftLimit*demandEcar(t,r);


genEqu(eq_plant_capacity,(t,rp(r,p)),(t,r,p))..
$ifthene.scaling NOT %SCALING% == 1
    genVar(POWER,(t,rp)) =l= (plant_cap(t,rp) + genVar(PLANT_ADD_CAP,(rp))*%RESOLUTION%) * avail(t,rp) ;
$else.scaling
    [genVar(POWER,(t,rp)) - (plant_cap(t,rp) + genVar(PLANT_ADD_CAP,(rp))*%RESOLUTION%) * avail(t,rp)] / [1$(avail(t,rp)<1e-6) + ((%RESOLUTION%)*avail(t,rp)*10)$(avail(t,rp)>=1e-6)] =l= 0;
$endif.scaling

genEqu(eq_total_plant_capacity,(rp(r,p)),(r,p))..
    sum(t, genVar(POWER,(t,rp))) + linkEqEpsDummy(eq_total_plant_capacity,(r,p)) =l= total_plant_cap(rp);

genEqu(eq_storage_balance,(tt,rs(r,s)),(tt,r,s))..
    genVar(STORAGE_LEVEL,(tt,rs))
  + linkEqEpsDummy(eq_storage_balance,(tt,r,s)) =e= genVar(STORAGE_LEVEL,(tt--1,rs))$t2(tt--1) * storage_efficiency(rs)
                                                  + genVar(STORAGE_INFLOW,(tt,rs))    * storage_efficiency_in(rs)
                                                  - genVar(STORAGE_OUTFLOW,(tt,rs))   / storage_efficiency_out(rs);

genEqu(eq_storage_capacity,(t,rs(r,s)),(t,r,s))..
    genVar(STORAGE_LEVEL,(t,rs)) =l= storage_cap(rs) + genVar(STORAGE_ADD_CAP,(rs));

genEqu(eq_emission_cap,(e),(e))..
$ifthene.scaling NOT %SCALING% == 1
    sum((rp(r,p),t), genVar(POWER,(t,rp)) * plant_emission(rp,e)) + linkEqEpsDummy(eq_emission_cap,(e)) =l= total_emission_cap(e);
$else.scaling
    [sum((rp(r,p),t), genVar(POWER,(t,rp)) * plant_emission(rp,e)) - total_emission_cap(e)]
    / [10*smin((rp(r,p),t)$plant_emission(rp,e), plant_emission(rp,e))]
    + linkEqEpsDummy(eq_emission_cap,(e))
    =l= 0;
$endif.scaling

genEqu(eq_link_capacity,(t,net),(t,net))..
    genVar(FLOW,(t,net)) + linkEqEpsDummy(eq_link_capacity,(t,net)) =l= link_cap(t,net) + genVar(LINK_ADD_CAP,(net)) * %RESOLUTION%;

model simple / all /;

$ifthene.noslack %NOSLACK%==1
  SLACK.fx(tt,rr)  = 0;
  simple.holdfixed = 1;
$endif.noslack


$iftheni.method %METHOD%==PIPS
   option lp=convertd;
   simple.optfile=1;
   simple.solvelink=2;
   simple.limrow=0;
   simple.limcol=0;
   simple.solprint=0;
*  annotation include file
$onechoV > %gams.scrdir%annotate.gms
   // ################ START ANNOTATION ########################################
   //Annotate variables
   DUMMY.stage                       = 1; DUMMY.lo = 0; Dummy.up = 1e-3;
   POWER.stage(t,rp(r,p))            = tr_stage(t,r);
   ECAR_BALANCE.stage(t,r,shiftclass)= tr_stage(t,r);
   ECAR_DELAY.stage(t,r,shiftclass)  = tr_stage(t,r);
   FLOW.stage(t,net(r1,r2))          = trr_stage(t,r1,r2);
   LINK_ADD_CAP.stage(net(r1,r2))$(card(tb)>=2) = 1;
   LINK_ADD_CAP.stage(net(r1,r2))$(card(tb)<=1) = smin(rtbmap(r1,ttmin,b(bb)), ord(bb))$[sum((rtbmap(r1,ttmin,b),rtbmap2(r2,ttmin,b)), 1)] + 1;
   SLACK.stage(t,r)                  = tr_stage(t,r);
   STORAGE_LEVEL.stage(t,rs(r,s))    = tr_stage(t,r);
   STORAGE_INFLOW.stage(t,rs(r,s))   = tr_stage(t,r);
   STORAGE_OUTFLOW.stage(t,rs(r,s))  = tr_stage(t,r);
   if(card(tb)<=1,
     PLANT_ADD_CAP.stage(rp(r,p))    = smin(rtbmap(r,tt,b(bb)), ord(bb)) + 1; // if no time decomposition put in region block
     STORAGE_ADD_CAP.stage(rs(r,s))  = smin(rtbmap(r,tt,b(bb)), ord(bb)) + 1; // if no time decomposition put in region block
   else
     PLANT_ADD_CAP.stage(rp(r,p))    =  1; // if time decomposition put in linking block
     STORAGE_ADD_CAP.stage(rs(r,s))  =  1; // if time decomposition put in linking block
   );
   BOBJ.stage(bb)                    = ord(bb) + 1;

   gen_POWER(t,rp(r,p))              = sum(trbmap(t,r,b),1);
   gen_ECAR_BALANCE(t,r,shiftclass)  = sum(trbmap(t,r,b),1);
   gen_ECAR_DELAY(t,r,shiftclass)    = sum(trbmap(t,r,b),1);
   gen_FLOW(t,net(r1,r2))            = not samerb(r1,r2) or sum(trbmap(t,r1,b),1);
   gen_LINK_ADD_CAP(net(r1,r2))      = yes;
   gen_SLACK(t,r)                    = yes;
   gen_STORAGE_LEVEL(t,rs(r,s))      = sum(trbmap(t,r,b),1);
   gen_STORAGE_INFLOW(t,rs(r,s))     = sum(trbmap(t,r,b),1);
   gen_STORAGE_OUTFLOW(t,rs(r,s))    = sum(trbmap(t,r,b),1);
   gen_PLANT_ADD_CAP(rp(r,p))        = yes;
   gen_STORAGE_ADD_CAP(rs(r,s))      = yes;
   gen_ROBJ(r)                       = yes;
   gen_BOBJ(b)                       = yes;
$  ifthene %rblock%>=2
     lvar_FLOW(tt,net(rr1,rr2))      = trr_stage(tt,net) = 1;
$  else
     lvar_FLOW(tt,net(rr1,rr2))      = no;
$  endif
   lvar_LINK_ADD_CAP(net(rr1,rr2))   = not samerb(rr1,rr2) or card(tb)>=2;
   lvar_PLANT_ADD_CAP(rp(rr,p))      = card(tb)>=2;
   lvar_STORAGE_ADD_CAP(rs(rr,s))    = card(tb)>=2;

   //annotate equations
   eq_bobj.stage(b(bb))               = ord(bb) + 1;
   eq_power_balance.stage(t,r)        = tr_stage(t,r);

   eq_ecar_balance.stage(t,r,shiftclass)$(not t_shift_linking(t,shiftclass))  = tr_stage(t,r);
   eq_ecar_balance.stage(tt,r,shiftclass)$t_shift_linking(tt,shiftclass) = (card(bb)+2);

   eq_ecar_shift_limit.stage(t,r)     = tr_stage(t,r);

   eq_plant_capacity.stage(t,rp(r,p)) = tr_stage(t,r);
   if(card(tb)<=1,
     eq_total_plant_capacity.stage(rp(r,p)) = rt_stage(r,ttmin); // no time decomposition --> region block
   else
     eq_total_plant_capacity.stage(rp(r,p)) = (card(bb)+2);     // time decomposition    --> coupling constraint
   );

   eq_storage_balance.stage(t,rs(r,s))$(not t_strg_linking(t)) = tr_stage(t,r); //non-linking
   eq_storage_balance.stage(t_strg_linking,rs(r,s)) = card(bb) + 2;             //linking

   eq_storage_capacity.stage(t,rs(r,s)) = tr_stage(t,r);
   eq_emission_cap.stage(e)               = card(bb) + 2;       //always linking
   eq_link_capacity.stage(t,net(r1,r2))   = trr_stage(t,r1,r2);
   eq_obj.stage = card(bb) + 2;

   gen_eq_bobj(b)                        = yes;
   gen_eq_power_balance(t,r)             = sum(trbmap(t,r,b), 1);

   gen_eq_ecar_balance(t,r,shiftclass)$(not t_shift_linking(t,shiftclass))  =  sum(trbmap(t,r,b), 1);
   gen_eq_ecar_balance(tt,r,shiftclass)$t_shift_linking(tt,shiftclass) = yes;

   gen_eq_ecar_shift_limit(t,r)          = sum(trbmap(t,r,b), 1);


   gen_eq_plant_capacity(t,rp(r,p))      = sum(trbmap(t,r,b), 1);
   gen_eq_total_plant_capacity(rp(r,p))  = sum(trbmap(t,r,b), 1) or card(tb)>=2;
   gen_eq_storage_balance(t,rs(r,s))     = sum(trbmap(t,r,b), 1); //non-linking constraints
   gen_eq_storage_balance(t_strg_linking,rs(r,s)) = yes;          //add the linking strg balance constraints
   gen_eq_storage_capacity(t,rs(r,s))    = sum(trbmap(t,r,b), 1);
   gen_eq_emission_cap(e)                = yes;
   gen_eq_link_capacity(t,net(r1,r2))    = sum((trbmap(t,r1,b),samerb(r1,r2)), 1) or (not samerb(r1,r2) and %BLOCK%=0);
   // ################ END ANNOTATION ########################################
$offecho
$onecho > %gams.scrdir%clearDomainSets.gms
*variable domain sets
option clear=gen_POWER, clear=gen_ECAR_BALANCE, clear=gen_ECAR_DELAY, clear=gen_FLOW, clear=gen_LINK_ADD_CAP, clear=gen_SLACK, clear=gen_STORAGE_LEVEL, clear=gen_STORAGE_INFLOW;
option clear=gen_STORAGE_OUTFLOW, clear=gen_PLANT_ADD_CAP, clear=gen_STORAGE_ADD_CAP, clear=gen_ROBJ, clear=gen_BOBJ;
option clear=lvar_FLOW, clear=lvar_LINK_ADD_CAP, clear=lvar_PLANT_ADD_CAP, clear=lvar_STORAGE_ADD_CAP;
*equation domain sets
option clear=gen_eq_bobj, clear=gen_eq_power_balance, clear=gen_eq_ecar_balance, clear=gen_eq_ecar_shift_limit, clear=gen_eq_plant_capacity, clear=gen_eq_total_plant_capacity, clear=gen_eq_storage_balance;
option clear=gen_eq_storage_balance, clear=gen_eq_storage_capacity, clear=gen_eq_emission_cap, clear=gen_eq_link_capacity;
$offecho
$  ifthene.blk %BLOCK%==-2
*    Generate single large gdx file
     t(tt) = yes;
     r(rr) = yes;
     b(bb) = yes;
$    include %gams.scrdir%annotate.gms
     put$(not %SUPPRESSDM%) fopt 'dictmap allblocks%GDXSUFFIX%_DM.gdx' /;
     putclose fopt 'jacobian allblocks%GDXSUFFIX%.gdx'
*	 execute_loadpoint 'allblocks%GDXSUFFIX%_sol'; option lp=default;
     solve simple min OBJ use lp;
*	 abort.noerror 'stop';
     execute 'mv -f %gams.scrdir%gamsdict.dat allblocks%GDXSUFFIX%_dict.gdx';

$    ifthene.runpips %RUNPIPS%==1     
$      eval bmaxp1 %bmax%+1            
       execute 'gmschk -t -X -g "%gams.sysdir%" %bmaxp1% allblocks%GDXSUFFIX%.gdx'
$      iftheni.computername %system.computername%==anton.gams.com
         execute 'mpirun -n %bmax% gmspips %bmaxp1% allblocks%GDXSUFFIX% "%gams.sysdir%" printsol'
$      else.computername
         execute 'srun -n %bmax% gmspips %bmaxp1% allblocks%GDXSUFFIX% "%gams.sysdir%" printsol'
$      endif.computername	   
       execute_loadpoint 'allblocks%GDXSUFFIX%_sol';
	   display '### OBJECTIVE FUNCTION VALUE:', OBJ.l;	   
$    endif.runpips

$  elseife.blk %BLOCK%==-1
*    Generate all (small) gdx files in a row
     equStage1 = card(bf) + 1; //enable generation of linking constraints for all block files
     loop(bf,
       actbf(bf) = yes;

       b(bb) = ord(bb) = ord(bf)-1;
       t(tt) = sum(trbmap(tt,rr,b),1) or (ord(bf)=1 and card(rb)>=2);
       r(rr) = yes; //for now we do not slice for regional decompoition
$      include %gams.scrdir%clearDomainSets.gms
$      include %gams.scrdir%annotate.gms
       put$(not %SUPPRESSDM%) fopt 'dictmap block%GDXSUFFIX%_DM' (ord(bf)-1):0:0 '.gdx' /;
       putclose fopt 'jacobian block%GDXSUFFIX%' (ord(bf)-1):0:0 '.gdx'
       solve simple min OBJ use lp;
*      put_utility 'exec' / 'mv -f %gams.scrdir%gamsdict.dat gamsdict' (ord(bf)-1):0:0 '.gdx';
       actbf(bf) = no;
     );

$  elseife.blk %BLOCK%>=0
*    Generate one particular (small) gdx file
     equStage1 = card(bf) + 1; //enable generation of linking constraints for all block files
     actbf(bf) = ord(bf) = %BLOCK%+1;

     b(bb) = ord(bb) = %BLOCK%; //define active block
     t(tt) = sum(trbmap(tt,rr,b),1) or (%BLOCK%=0 and card(rb)>=2); //define active time steps
     r(rr) = yes; //for now we do not slice for regional decompoition

$    include %gams.scrdir%annotate.gms
     put$(not %SUPPRESSDM%) fopt 'dictmap block%GDXSUFFIX%_DM%BLOCK%.gdx' /;
     putclose fopt 'jacobian block%GDXSUFFIX%%BLOCK%.gdx';
     solve simple min OBJ use lp;
*     execute 'mv -f %gams.scrdir%gamsdict.dat gamsdict%BLOCK%.gdx';
     actbf(bf) = no;

$  else.blk
     abort 'invalid setting of --BLOCK
$  endif.blk

$ELSE.method
*  Standard LP
   t(tt) = yes;
   r(rr) = yes;

   option clear=eq_bobj, clear=BOBJ, clear=POWER, //clear=SLACK,
          clear=eq_power_balance, clear=eq_ecar_balance,
          clear=FLOW, clear=eq_plant_capacity, clear=eq_total_plant_capacity, clear=eq_storage_balance,
          clear=STORAGE_LEVEL, clear=eq_storage_capacity, clear=eq_emission_cap, clear=eq_link_capacity,
          clear=eq_obj;
   option clear=lvar_FLOW, clear=lvar_LINK_ADD_CAP, clear=lvar_PLANT_ADD_CAP, clear=lvar_STORAGE_ADD_CAP;

   gen_POWER(t,rp)                      = yes;
   gen_ECAR_BALANCE(t,r,shiftclass)     = yes;
   gen_ECAR_DELAY(t,r,shiftclass)       = yes;
   gen_FLOW(t,net(r1,r2))               = yes;
   gen_LINK_ADD_CAP(net(r1,r2))         = yes;
   gen_SLACK(t,r)                       = yes;
   gen_STORAGE_LEVEL(t,rs(r,s))         = yes;
   gen_STORAGE_INFLOW(t,rs(r,s))        = yes;
   gen_STORAGE_OUTFLOW(t,rs(r,s))       = yes;
   gen_PLANT_ADD_CAP(rp(r,p))           = yes;
   gen_STORAGE_ADD_CAP(rs(r,s))         = yes;
   gen_ROBJ(r)                          = yes;
   gen_BOBJ(b)                          = yes;
   gen_eq_bobj(b)                       = yes;
   gen_eq_power_balance(t,r)            = yes;
   gen_eq_ecar_balance(t,r,shiftclass)  = yes;
   gen_eq_ecar_shift_limit(t,r)         = yes;
   gen_eq_plant_capacity(t,rp(r,p))     = yes;
   gen_eq_total_plant_capacity(rp(r,p)) = yes;
   gen_eq_storage_balance(t,rs(r,s))    = yes;
   gen_eq_storage_capacity(t,rs(r,s))   = yes;
   gen_eq_emission_cap(e)               = yes;
   gen_eq_link_capacity(t,net(r1,r2))   = yes;

   simple.optfile = 1;
$onecho > cplex.opt
lpmethod 4
solutiontype 2
*preind 0
$offecho

*ECAR_BALANCE.fx(t,r,shiftclass)=0;
*ECAR_DELAY.fx(t,r,shiftclass)=0;
execute_loadpoint '%sol%'; simple.limrow=200000; simple.limcol=20000;
   solve simple min OBJ use lp;

$ENDIF.method


