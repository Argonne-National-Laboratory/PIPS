$ eolcom //

$if  not set FROM                 $set FROM              0
$if  not set TO                   $set TO                1
$if  not set RESOLUTION           $set RESOLUTION        1
$if  not set NBREGIONS            $set NBREGIONS         4
$if  not set METHOD               $set METHOD            standard_lp
$if  not set EMPMETHOD            $set EMPMETHOD         DE

$if  not set LOADFROMXLS          $ set LOADFROMXLS      0
$if  not set XLSID                $ set XLSID            standard

$if  not set SCALING              $set SCALING           0
$if  not set NOSLACK              $set NOSLACK           0


$ife %FROM%>%TO%         $abort 'FROM > TO'
$ife %RESOLUTION%<0      $abort 'Negative RESOLUTION forbidden'

$eval NBTIMESTEPS ceil((%TO%-%FROM%)*365*24/%RESOLUTION%)
$setnames "%gams.input%" SIMPLEDIR fn fe

* if there is no input data for the current settings available call data generator
$ifthen.xls %LOADFROMXLS%==0
$  include "%SIMPLEDIR%simple_data_gen.gms" 
$else.xls
$  call gams "%SIMPLEDIR%simple_data_excel.gms" --XLSID=%XLSID%
$  if errorlevel 1 $abort 'problems with creating REMix data from Excel'
$endif.xls


*   basic sets
Set rr                          'regions                                       '
    rrUel
    p                           'plants                                        '
    s                           'storages                                      '
    tt                          'time steps                                    ' / t1*t%NBTIMESTEPS% /
    ttX                         'time steps in input data                      '
    e                           'emissions                                     '
    type                        'plant type                                    '
;

Alias(rr,rr1,rr2,rr3);

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
Alias(r,r1,r2);
$iftheni.method %METHOD%==spExplicitDE
$ if not set NBSCEN $set NBSCEN 5
  set scen / scen1*scen%NBSCEN% /
      sc(scen);
  parameter prob(scen);
$ set SCEN scen,
$ set SCENS sc,
$else.method
$ set SCEN
$ set SCENS
$endif.method

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
    type_mult(%SCEN%type)        'generation cost multiplier                   '
    availX(ttX,rr,p)
    avail(tt,rr,p)
    avail2(rr,p,tt)
;

scalar
    tmpStart                     'helper storing time step start times         '
    tmpEnd                       'helper storing time step end times           '
;

* read data from gdx
$ifthen.xls %LOADFROMXLS%==1
$ gdxin "%SIMPLEDIR%REMixData_%XLSID%.gdx"
$ load ttX
$ load rrUel=rr
$ iftheni.method %METHOD%==pipsX
$  ifthene.blk %BLOCK%>0
$   onecho > genRegByBlock.gms
    set rr; alias (rr,rr1,rr2); set net(rr1,rr2);
$   ifthen %LOADFROMXLS%==0
$    gdxin "%SIMPLEDIR%simple_data_%NBREGIONS%.gdx"
$   else
$    gdxin "%SIMPLEDIR%REMixData_%XLSID%.gdx"
$   endif
$   load rr net
    set rrBlock(rr); singleton set r(rr); r(rr) = ord(rr)=%BLOCK%;
    rrBlock(rr) = net(rr,r) or net(r,rr); rrBlock(r) = yes;
    execute_unload '%SIMPLEDIR%genRegByBlock', rrBlock;
$   offecho
$   call gams genRegByBlock lo=%gams.lo%
$   if errorlevel 1 $abort problems with genRegByBlock
$   gdxin genRegByBlock
$   load rr=rrBlock
$   ifthen.xls2 %LOADFROMXLS%==0
$    gdxin "%SIMPLEDIR%simple_data_%NBREGIONS%.gdx"
$   else.xls2
$    gdxin REMixData_%XLSID%.gdx
$   endif.xls2
$  else.blk
$   load rr
$  endif.blk
$ else.method
$  load rr
$ endif.method
$ load p s e pr rp type ptype re sr rs net
$ load plant_capX cost_unserved_demandX link_capX link_efficiencyX demandX
$ load yearly_plant_cap plant_emission cost_power_generation yearly_emission_cap
$ load cost_emission storage_cap storage_efficiency storage_efficiency_in
$ load storage_efficiency_out storage_max_in storage_max_out
$ load plant_max_add_cap storage_max_add_cap link_max_add_cap cost_plant_add
$ load cost_storage_add cost_link_add availX
$ gdxin
$endif.xls

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

Scalar genBlock indicator for generating a block in case of method=PIPS /0/;


Positive variables
    POWER(%SCEN%tt,rr,p)           'power production                   [GWh]'
$IFTHENI.method %METHOD%==lagrange_relaxation
    FLOW(tt,rr,rr1,rr2)            'power flow                         [GWh]'
    LINK_ADD_CAP(rr,rr1,rr2)       'arc capacity expansion             [GWh]'
$ELSE.method
    FLOW(%SCEN%tt,rr1,rr2)         'power flow                         [GWh]'
    LINK_ADD_CAP(rr1,rr2)          'arc capacity expansion             [GWh]'
$ENDIF.method
    SLACK(%SCEN%tt,rr)             'uncovered demand                   [GWh]'
    STORAGE_LEVEL(%SCEN%tt,rr,s)   'storage level                      [GWh]'
    STORAGE_INFLOW(%SCEN%tt,rr,s)  'power entering storage             [GWh]'
    STORAGE_OUTFLOW(%SCEN%tt,rr,s) 'power taken from storage           [GWh]'
    PLANT_ADD_CAP(rr,p)            'plant capacity expansion            [GW]'
    STORAGE_ADD_CAP(rr,s)          'storage capacity expansion         [GWh]'
    EMISSION_SPLIT(%SCEN%rr,e)     'emission allowance split      [fraction]'
Variable
    EMISSION_COST(%SCEN%rr,e)      'emission cost                     [MEUR]'
    ROBJ(%SCEN%rr)                 'total region cost                 [MEUR]'
    OBJ                            'total                             [MEUR]'
;

* bounds
$IFTHENI.method %METHOD%==lagrange_relaxation
netx(rr1,net(rr1,rr2))          = yes;
netx(rr2,net(rr1,rr2))          = yes;
LINK_ADD_CAP.up(netx(rr,net))   = link_max_add_cap(net)  ;
$ELSE.method
LINK_ADD_CAP.up(net)            = link_max_add_cap(net)  ;
$ENDIF.method
PLANT_ADD_CAP.up(rp)            = plant_max_add_cap(rp)  ;
STORAGE_ADD_CAP.up(rs)          = storage_max_add_cap(rs);
* The bounds for STORAGE_INFLOW and STORAGE_OUTFLOW are set in the include file spexplicitde.gms
$iftheni.method not %METHOD%==spExplicitDE
  STORAGE_INFLOW.up(tt,rs)  = storage_max_in(rs)     ;
  STORAGE_OUTFLOW.up(tt,rs) = storage_max_out(rs)    ;
$endif.method

equations
    eq_robj(%SCEN%rr)                     'total cost in region                 '
    eq_power_balance(%SCEN%tt,rr)         'power balance                        '
    eq_plant_capacity(%SCEN%tt,rr,p)      'respect plant capacity               '
    eq_total_plant_capacity(%SCEN%rr,p)   'respect total plant capacity         '
    eq_storage_balance(%SCEN%tt,rr,s)     'storage balance                      '
    eq_storage_capacity(%SCEN%tt,rr,s)    'respect storage capacity             '
    eq_emission_region(%SCEN%rr,e)        'calculate regional emissions         '
    eq_emission_cost(%SCEN%rr,e)          'calculate regional emission costs    '
    eq_emission_cap(%SCEN%e)              'respect emission cap                 '
$IFTHENI.method %METHOD%==lagrange_relaxation
    eq_link_capacity(tt,rr,rr1,rr2)       'respect link capacity                '
    eq_same_add_cap(rr1,rr2)              'same capacity extension on copies    '
    eq_same_flow(tt,rr1,rr2)              'same flow on copies of flow variables'
$ELSE.method
    eq_link_capacity(%SCEN%tt,rr1,rr2)    'respect link capacity                '
$ENDIF.method
    eq_obj                                'total cost                           '
;

eq_obj..
    OBJ =e=
$IFTHENI.method %METHOD%==spExplicitDE
            sum(sc, prob(sc)*sum(r, ROBJ(sc,r)))
          + sum(rp(r,p),     PLANT_ADD_CAP(rp)              * EPS)
          + sum(rs(r,s),     STORAGE_ADD_CAP(rs)            * EPS)
$ELSEIFI.method not %METHOD%==distributedsimplebendersstage
            sum(r,           ROBJ(%SCENS%r))
$ENDIF.method
$IFTHENI.method %METHOD%==pips
          + sum((t,net(rr1,rr2)),  eps*FLOW(t,net))
          + sum((rr1,e),  eps*EMISSION_SPLIT(rr1,e))
$ENDIF.method
$IFTHENI.method %METHOD%==lagrange_relaxation
          + sum(netx(r,net(rr1,rr2)),  LINK_ADD_CAP(netx)   * cost_link_add(net))
$ELSEIFI.method %METHOD%==pipsx
          + sum(net(rr1,rr2),  LINK_ADD_CAP(net)   * cost_link_add(net))$(not genBlock)
          + sum(net(r,rr2),  LINK_ADD_CAP(net)   * cost_link_add(net))$genBlock
          + sum(net(rr2,r),  LINK_ADD_CAP(net)   * cost_link_add(net))$genBlock
$ELSE.method
          + sum(net(rr1,rr2),  LINK_ADD_CAP(net)   * cost_link_add(net))
$ENDIF.method
          ;
eq_robj(%SCENS%r)..
$IFTHENE.scaling NOT %SCALING% == 1
$  IFTHENI.method not %METHOD%==distributedsimplebendersstage
    ROBJ(%SCENS%r)
$  else.method
    OBJ
$  endif.method
            =e= sum((t,ptype(rp(r,p),type)), POWER(%SCENS%t,rp) * cost_power_generation(rp) * type_mult(%SCENS%type))
              + sum(t,           SLACK(%SCENS%t,r)              * cost_unserved_demand(t))
              + sum(rp(r,p),     PLANT_ADD_CAP(rp)              * cost_plant_add(rp))
              + sum(rs(r,s),     STORAGE_ADD_CAP(rs)            * cost_storage_add(rs))
              + sum(e,           EMISSION_COST(%SCENS%r,e));
$ELSE.scaling
$  IFTHENI.method not %METHOD%==distributedsimplebendersstage
    [ROBJ(%SCENS%r)
$  else.method
    [OBJ
$  endif.method
   - ( sum((t,ptype(rp(r,p),type)), POWER(%SCENS%t,rp) * cost_power_generation(rp) * type_mult(%SCENS%type))
     + sum(t,           SLACK(%SCENS%t,r)              * cost_unserved_demand(t))
     + sum(rp(r,p),     PLANT_ADD_CAP(rp)              * cost_plant_add(rp))
     + sum(rs(r,s),     STORAGE_ADD_CAP(rs)            * cost_storage_add(rs))
     + sum(e,           EMISSION_COST(%SCENS%r,e))) ]
   / [10*smin((t,ptype(rp(r,p),type))$cost_power_generation(rp), cost_power_generation(rp) * type_mult(%SCENS%type))]
   =e= 0;

$ENDIF.scaling

eq_power_balance(%SCENS%t,r)..
        sum(rp(r,p),    POWER(%SCENS%t,rp))
$IFTHENI.method %METHOD%==lagrange_relaxation
      + sum(net(rr2,r), FLOW(t,r,rr2,r) * link_efficiency(t,net))
      - sum(net(r,rr2), FLOW(t,r,r,rr2))
$ELSE.method
      + sum(net(rr2,r), FLOW(%SCENS%t,net) * link_efficiency(t,net))
      - sum(net(r,rr2), FLOW(%SCENS%t,net))
$ENDIF.method
      + sum(rs(r,s),    STORAGE_OUTFLOW(%SCENS%t,rs) - STORAGE_INFLOW(%SCENS%t,rs))
      + SLACK(%SCENS%t,r)
    =g= demand(t,r);

eq_plant_capacity(%SCENS%t,rp(r,p))..
$ifthene.scaling NOT %SCALING% == 1
    POWER(%SCENS%t,rp) =l= (plant_cap(t,rp) + PLANT_ADD_CAP(rp)*%RESOLUTION%) * avail(t,rp) ;
$else.scaling
    [POWER(%SCENS%t,rp) - (plant_cap(t,rp) + PLANT_ADD_CAP(rp)*%RESOLUTION%) * avail(t,rp)] / [1$(avail(t,rp)<1e-6) + ((%RESOLUTION%)*avail(t,rp)*10)$(avail(t,rp)>=1e-6)] =l= 0;
$endif.scaling


eq_total_plant_capacity(%SCENS%rp(r,p))..
    sum(t, POWER(%SCENS%t,rp)) =l= total_plant_cap(rp);

eq_storage_balance(%SCENS%t(tt),rs(r,s))..
    STORAGE_LEVEL(%SCENS%t,rs) =e= STORAGE_LEVEL(%SCENS%tt--1,rs) * storage_efficiency(rs)
                         + STORAGE_INFLOW(%SCENS%t,rs)    * storage_efficiency_in(rs)
                         - STORAGE_OUTFLOW(%SCENS%t,rs)   / storage_efficiency_out(rs) ;

eq_storage_capacity(%SCENS%t,rs(r,s))..
    STORAGE_LEVEL(%SCENS%t,rs) =l= storage_cap(rs) + STORAGE_ADD_CAP(rs);

eq_emission_region(%SCENS%r,e)..
$ifthene.scaling NOT %SCALING% == 1
    sum((rp(r,p),t), POWER(%SCENS%t,rp) * plant_emission(rp,e)) =l= total_emission_cap(e)*EMISSION_SPLIT(%SCENS%r,e);
$else.scaling
    [sum((rp(r,p),t), POWER(%SCENS%t,rp) * plant_emission(rp,e)) - total_emission_cap(e)*EMISSION_SPLIT(%SCENS%r,e)]
    / [10*smin((rp(r,p),t)$plant_emission(rp,e), plant_emission(rp,e))] =l= 0;
$endif.scaling

eq_emission_cost(%SCENS%r,e)..
$ifthene.scaling NOT %SCALING% == 1
    sum((rp(r,p),t), POWER(%SCENS%t,rp) * plant_emission(rp,e)) * cost_emission(e) =e= EMISSION_COST(%SCENS%r,e);
$else.scaling
      [sum((rp(r,p),t), POWER(%SCENS%t,rp) * plant_emission(rp,e)) * cost_emission(e) - EMISSION_COST(%SCENS%r,e)]
    / [10*smin((rp(r,p),t)$plant_emission(rp,e), plant_emission(rp,e)*cost_emission(e))]   =e= 0;
$endif.scaling

eq_emission_cap(%SCENS%e)$(not genBlock)..
    sum(rr, EMISSION_SPLIT(%SCENS%rr,e)) =l= 1;

$IFTHENI.method %METHOD%==lagrange_relaxation
eq_link_capacity(t,netx(r,net(rr1,rr2)))..
    FLOW(t,netx) =l= link_cap(t,net) + LINK_ADD_CAP(netx) * %RESOLUTION%;

eq_same_flow(t,net(rr1,rr2))..
    FLOW(t,rr2,rr1,rr2) =e= FLOW(t,rr1,rr1,rr2);

eq_same_add_cap(net(rr1,rr2))..
    LINK_ADD_CAP(rr2,rr1,rr2) =e= LINK_ADD_CAP(rr1,rr1,rr2);

$ELSE.method
eq_link_capacity(%SCENS%t,net)$(not genBlock)..
    FLOW(%SCENS%t,net) =l= link_cap(t,net) + LINK_ADD_CAP(net) * %RESOLUTION%;
$ENDIF.method

$IFTHENI %METHOD%_%EMPMETHOD%==stochasticEMP_LindoBenders
equation dummyForLindoBenders;
dummyForLindoBenders.. 1 =g= 0;
$endif

model simple / all /;
option limrow=0, limcol=0, solprint=silent;

$ifthene.noslack %NOSLACK%==1
$ IFTHENI.method %METHOD%==spExplicitDE
   sc(scen) = yes;
$ ENDIF.method
  SLACK.fx(%SCENS%tt,rr)  = 0;
  simple.holdfixed = 1;
$ IFTHENI.method %METHOD%==spExplicitDE
   sc(scen) = no;
$ ENDIF.method
$endif.noslack

$IFTHENI.method %METHOD%==rolling_horizon
   r(rr) = yes;
$  include "%SIMPLEDIR%rolling_horizon.gms"

$ELSEIFI.method %METHOD%==spExplicitDE
$  include "%SIMPLEDIR%spexplicitde.gms"

$ELSEIFI.method %METHOD%==benders
$  include "%SIMPLEDIR%benders.gms"

$ELSEIFI.method %METHOD%==bendersMPI
$  include "%SIMPLEDIR%bendersMPI.gms"

$ELSEIFI.method %METHOD%==bendersMPImi
$  include "%SIMPLEDIR%bendersMPImi.gms"

$ELSEIFI.method %METHOD%==simplebenders
   file fbaf / simple.baf /; put fbaf;
$  onput
   masterVar:  FLOW LINK_ADD_CAP EMISSION_SPLIT
   masterEqu:  eq_emission_cap eq_link_capacity
   subVar:     ROBJ#1 POWER#2 SLACK#2 STORAGE_LEVEL#2 STORAGE_INFLOW#2
               STORAGE_OUTFLOW#2 PLANT_ADD_CAP#1 STORAGE_ADD_CAP#1
               EMISSION_COST#1
   subEqu:     eq_robj#1 eq_power_balance#2 eq_plant_capacity#2
               eq_total_plant_capacity#1 eq_storage_balance#2
               eq_storage_capacity#2 eq_emission_region#1 eq_emission_cost#1
   blockIndex:
$  offput
   loop(rr, put rr.tl:0 " ");
   putclose;
   t(tt) = yes;
   r(rr) = yes;
   FLOW.l(t,net)          = 0;
   LINK_ADD_CAP.l(net)    = 0;
   Scalar totalDemand; totalDemand = sum((rr2,t), demand(t,rr2));
   EMISSION_SPLIT.l(r,e)$totalDemand = sum(t,demand(t,r)) / totalDemand;
$ontext
In order to have simplebenders (and others) work one needs to create a file
simplecmpNT.txt with the folloing line and start the GAMS job with subsys=simplecmpNT.txt:

CONVERTD 102011 15 0001020304 1 1 2 LP
gmsgennt.cmd
gmsgennx.exe
cvddclib.dll cvd 1 0

DSBSTAGE 2001 15 0001020304 1 0 1 LP
noscript.cmd
C:\Users\Michael\Documents\share\fred\BEAM-ME\SIMPLE\SimpleBenders\Debug\DistributedSimpleBendersStage.exe

SB 1 5 0001020304 1 0 1 LP
noscript.cmd
C:\Users\Michael\Documents\share\fred\BEAM-ME\SIMPLE\SimpleBenders\Debug\SimpleBenders.exe

SBSTAGE 2001 5 0001020304 1 0 1 LP
noscript.cmd
C:\Users\Michael\Documents\share\fred\BEAM-ME\SIMPLE\SimpleBenders\Debug\SimpleBendersStage.exe

DEFAULTS
LP CONVERTD
$offtext
   option lp=sb, solvelink=2;
   solve simple min OBJ use lp;

$ELSEIFI.method %METHOD%==simplebendersstage
   t(tt) = yes;
   r(rr) = yes;
   FLOW.l(t,net)          = 0;
   LINK_ADD_CAP.l(net)    = 0;
   Scalar totalDemand; totalDemand = sum((rr2,t), demand(t,rr2));
   EMISSION_SPLIT.l(r,e)$totalDemand = sum(t,demand(t,r)) / totalDemand;
* Master variables and equation
   FLOW.stage(t,net) = 1;
   LINK_ADD_CAP.stage(net) = 1;
   EMISSION_SPLIT.stage(rr,e) = 1;
   eq_emission_cap.stage(e) = 1;
   eq_link_capacity.stage(t,net) = 1;
* Block variables and equations
   ROBJ.stage(rr) = ord(rr)+1;
   POWER.stage(t,rp(rr,p)) = ord(rr)+1;
   SLACK.stage(t,rr) = ord(rr)+1;
   STORAGE_LEVEL.stage(t,rs(rr,s)) = ord(rr)+1;
   STORAGE_INFLOW.stage(t,rs(rr,s)) = ord(rr)+1;
   STORAGE_OUTFLOW.stage(t,rs(rr,s)) = ord(rr)+1;
   STORAGE_ADD_CAP.stage(rs(rr,s)) = ord(rr)+1;
   PLANT_ADD_CAP.stage(rp(rr,p)) = ord(rr)+1;
   EMISSION_COST.stage(rr,e) = ord(rr)+1;
   eq_robj.stage(rr) = ord(rr)+1;
   eq_power_balance.stage(t,rr) = ord(rr)+1;
   eq_plant_capacity.stage(t,rp(rr,p)) = ord(rr)+1;
   eq_total_plant_capacity.stage(rp(rr,p)) = ord(rr)+1;
   eq_storage_balance.stage(t,rs(rr,s)) = ord(rr)+1;
   eq_storage_capacity.stage(t,rs(rr,s)) = ord(rr)+1;
   eq_emission_region.stage(rr,e) = ord(rr)+1;
   eq_emission_cost.stage(rr,e) = ord(rr)+1;

$set DOCONVERT 1
$ifthen.convert set DOCONVERT
*$ontext
$  echo jacobian > convertd.opt
$onecho > subsys.txt
CONVERTD 102011 15 0001020304 1 1 2 LP
gmsgennt.cmd
gmsgennx.exe
cvddclib.dll cvd 1 0
$offecho

   simple.optfile=1;
   option lp=convertd;
*$offtext
$else.convert
   option lp=sbStage, solvelink=2;
$endif.convert
   solve simple min OBJ use lp;

$ELSEIFI.method %METHOD%==distributedsimplebendersstage
$if not set BLOCKMAX $eval BLOCKMAX card(rr)
set block /block1*block%BLOCKMAX%/;
$ifthen.block not set BLOCK
   t(tt) = yes;
   r(rr) = yes;
$else.block
$ontext
* Region and time decomposition
set blockinfo(block,rr,tt) /
   block1.(r1,r2).(t1*t168)
   block2.(r1,r2).(t169*t336)
   block3.(r3,r4).(t1*t168)
   block4.(r3,r4).(t169*t336)
/
$offtext
$ontext
* Time decomposition
set blockinfo(block,rr,tt) /
   block1.(#rr).(t1*t168)
   block2.(#rr).(t169*t336)
   block3.(#rr).(t337*t504)
   block4.(#rr).(t505*t672)
/
$offtext
*$ontext
* Region decomposition
set blockinfo(block,rr,tt) / (#block:#rr).#tt /;
*$offtext
set brt(rr,tt) dynamic slices for region and time;
   brt(rr,tt) = blockinfo('block%BLOCK%',rr,tt);
   option r<brt; option t<brt;
$endif.block
   FLOW.l(t,net)          = 0;
   LINK_ADD_CAP.l(net)    = 0;
   Scalar totalDemand; totalDemand = sum((rr2,t), demand(t,rr2));
   EMISSION_SPLIT.l(r,e)$totalDemand = sum(t,demand(t,r)) / totalDemand;
* Master variables and equation
   FLOW.stage(t,net) = 1;
   LINK_ADD_CAP.stage(net) = 1;
   EMISSION_SPLIT.stage(rr,e) = 1;
   eq_emission_cap.stage(e) = 1;
   eq_link_capacity.stage(t,net) = 1;
$ifthen.block set BLOCK
* Block variables and equations
   ROBJ.stage(r) = 2;
   POWER.stage(t,rp(r,p)) = 2;
   SLACK.stage(t,r) = 2;
   STORAGE_LEVEL.stage(t,rs(r,s)) = 2;
   STORAGE_INFLOW.stage(t,rs(r,s)) = 2;
   STORAGE_OUTFLOW.stage(t,rs(r,s)) = 2;
   STORAGE_ADD_CAP.stage(rs(r,s)) = 2;
   PLANT_ADD_CAP.stage(rp(r,p)) = 2;
   EMISSION_COST.stage(r,e) = 2;
   eq_robj.stage(r) = 2;
   eq_power_balance.stage(t,r) = 2;
   eq_plant_capacity.stage(t,rp(r,p)) = 2;
   eq_total_plant_capacity.stage(rp(r,p)) = 2;
   eq_storage_balance.stage(t,rs(r,s)) = 2;
   eq_storage_capacity.stage(t,rs(r,s)) = 2;
   eq_emission_region.stage(r,e) = 2;
   eq_emission_cost.stage(r,e) = 2;
$endif.block

   model simpleFirstStage /  eq_link_capacity, eq_emission_cap, eq_obj /;
   model simpleBlock /  simple - simpleFirstStage /;
$ifthen.block not set BLOCK
   option lp=dsbStage, solvelink=2;
* Pass on the number of blocks in integer5
   simpleFirstStage.integer5 = %BLOCKMAX%;

* Parameter file for GAMS calls from inside dsbStage
$onecho  > dsbStage.pf
input="%gams.input%"
subsys=simplecmpNT.txt
lo=4
--METHOD=%METHOD%
--TO=%TO% --FROM=%FROM% --RESOLUTION=%RESOLUTION% --NBREGIONS=%NBREGIONS%
$offecho

   solve simpleFirstStage min OBJ use lp;
   scalar blk; file fx;
   loop(rr,
      put_utility fx "gdxin" / "block" ord(rr):0:0 ".gdx";
      execute_loadpoint;
      robj.l(rr) = obj.l;
   )
   execute_loadpoint "master.gdx";
$else.block
   option lp=convertd, solvelink=2;
* Instruct to generate scratch files but that's it
   simpleBlock.integer1 = 234;
   solve simpleBlock min OBJ use lp;
$endif.block

$ELSEIFI.method %METHOD%==pips
   Scalar BLOCK0 /1/, BLOCKN;
   BLOCKN = BLOCK0 + card(rr);

   t(tt) = yes; r(rr) = no;
* Master variables and equation
   FLOW.stage(t,net) = BLOCK0;
   LINK_ADD_CAP.stage(net) = BLOCK0;
   EMISSION_SPLIT.stage(rr,e) = BLOCK0;
   eq_emission_cap.stage(e) = BLOCK0;
   eq_link_capacity.stage(t,net) = BLOCK0;

   file fopt / "%gams.optdir%convertd.opt" /;
$onecho > subsys.txt
CONVERTD 102011 15 0001020304 1 1 2 LP
gmsgennt.cmd
gmsgennx.exe
cvddclib.dll cvd 1 0
$offecho
   option lp=convertd;
   simple.optfile=1;
   simple.solvelink=2;
   simple.limrow=0;
   simple.limcol=0;
   simple.solprint=0;

$if not set BLOCK $set BLOCK -1
$ifthene.blk (%BLOCK%=0)or(%BLOCK%=-1)
* Block 0
   putclose fopt 'jacobian block0.gdx' / 'dictmap dmblock0.gdx';
   solve simple min OBJ use lp;
$endif.blk

   loop(sameas(rr3,rrUel)$((%BLOCK%<0) or (%BLOCK%=ord(rrUel))),
      option clear=r;
      r(rr3) = yes;
* Block variables and equations
      ROBJ.stage(rr3) = ord(rrUel)+BLOCK0;
      POWER.stage(t,rp(rr3,p)) = ord(rrUel)+BLOCK0;
      SLACK.stage(t,rr3) = ord(rrUel)+BLOCK0;
      STORAGE_LEVEL.stage(t,rs(rr3,s)) = ord(rrUel)+BLOCK0;
      STORAGE_INFLOW.stage(t,rs(rr3,s)) = ord(rrUel)+BLOCK0;
      STORAGE_OUTFLOW.stage(t,rs(rr3,s)) = ord(rrUel)+BLOCK0;
      STORAGE_ADD_CAP.stage(rs(rr3,s)) = ord(rrUel)+BLOCK0;
      PLANT_ADD_CAP.stage(rp(rr3,p)) = ord(rrUel)+BLOCK0;
      EMISSION_COST.stage(rr3,e) = ord(rrUel)+BLOCK0;
      eq_robj.stage(rr3) = ord(rrUel)+BLOCK0;
      eq_power_balance.stage(t,rr3) = ord(rrUel)+BLOCK0;
      eq_plant_capacity.stage(t,rp(rr3,p)) = ord(rrUel)+BLOCK0;
      eq_total_plant_capacity.stage(rp(rr3,p)) = ord(rrUel)+BLOCK0;
      eq_storage_balance.stage(t,rs(rr3,s)) = ord(rrUel)+BLOCK0;
      eq_storage_capacity.stage(t,rs(rr3,s)) = ord(rrUel)+BLOCK0;
      eq_emission_region.stage(rr3,e) = ord(rrUel)+BLOCK0;
      eq_emission_cost.stage(rr3,e) = ord(rrUel)+BLOCK0;
$ifthene.blk %BLOCK%>-2
      putclose fopt 'jacobian block' ord(rrUel):0:0 '.gdx' / 'dictmap dmblock' ord(rrUel):0:0 '.gdx';
      genBlock = 1;
      solve simple min OBJ use lp;
      genBlock = 0;
$endif.blk
   );
$ifthene.blk %BLOCK%=-2
      putclose fopt 'jacobian allblocks.gdx' / 'dictmap dmallblocks.gdx';
      r(rr) = yes;
      solve simple min OBJ use lp;
$endif.blk
$ELSEIFI.method %METHOD%==lagrange_relaxation
$  if not set lgrhpp $set lgrhpp 1
** solve standard LP
*   t(tt) = yes; r(rr) = yes; solve simple min OBJ use lp; t(tt) = no; r(rr) = no;
** use RH to get good first lambda
$ifthen.findLambda %lgrhpp%==1
   r(rr) = yes;
$  include "%SIMPLEDIR%rolling_horizon.gms"
   r(rr) = no;
$endif.findLambda
** include lagrange relaxation
   t(tt) = yes;
$  include "%SIMPLEDIR%lg.gms"

$ELSEIFI.method %METHOD%==stochasticEMP
$  include "%SIMPLEDIR%stochasticEMP.gms"

$ELSEIFI.method %METHOD%==SPBENDERSSEQ
   t(tt) = yes;
   r(rr) = yes;
$  include "%SIMPLEDIR%spBendersSeq.gms"

$ELSEIFI.method %METHOD%==SPBENDERSASYNC
   t(tt) = yes;
   r(rr) = yes;
$  include "%SIMPLEDIR%spBendersAsync.gms"

$ELSEIFI.method %METHOD%==SPBENDERSMPI
   t(tt) = yes;
   r(rr) = yes;
$  include "%SIMPLEDIR%spBendersMPI.gms"

$ELSE.method
*  Standard LP
   t(tt) = yes;
   r(rr) = yes;
   if (card(rr)>10 or card(tt)>50,
      option limrow=0, limcol=0;
   );
   simple.optfile = 1;
$  onecho > cplex.opt
   lpmethod     4
   solutiontype 2
$  offecho
   solve simple min OBJ use lp;
$ENDIF.method


