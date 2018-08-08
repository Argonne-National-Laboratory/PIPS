* Rolling Horizon Implementation by S. Schreck
$if  not set NBINTERVALSRH        $ set NBINTERVALSRH     4
$if  not set OVERLAPSIZERH        $ set OVERLAPSIZERH     0.1

option solvelink = 5, limrow = 0, limcol = 0, solprint = silent;

* solve coarse problem version first
$ontext
scalar coarse_resolution;
coarse_resolution = %NBTIMESTEPS% / %NBINTERVALSRH%;
file fx; put fx;
put_utility 'exec' / 'gams simple.gms lo=3 o=coarse.lst gdx=coarse.gdx --FROM=%FROM% --TO=%TO% --RESOLUTION=' coarse_resolution;
set tt_coarse(*), t_coarse(*);
parameter coarse_storage_level1(tt,rr,s)
          coarse_storage_level2(t_coarse,rr,s);
execute_load 'coarse.gdx', tt_coarse=tt, coarse_storage_level1=STORAGE_LEVEL.l;
display tt_coarse, coarse_storage_level1;
$offtext

***generating subsets for rolling-horizon***
set window  / w1*w%NBINTERVALSRH% /
set window_t(window,tt);
set window_t_fix(window,tt);
*alias (timeModelInterval,t)
scalars  first_hour, last_hour
         interval_size, overlap_size,
         interval_start, interval_end, interval_fix_end;


*first_hour = 1;
interval_fix_end = 0;
last_hour        = card(tt);
STORAGE_LEVEL.fx('t%NBTIMESTEPS%',rs) = 0.0 * storage_cap(rs);

overlap_size  = round(%OVERLAPSIZERH% * (%NBTIMESTEPS% / %NBINTERVALSRH%));
interval_size = round((%NBTIMESTEPS% - overlap_size) / %NBINTERVALSRH% + overlap_size);

loop(window,
  interval_start   = interval_fix_end + 1;
  interval_end     = interval_start + interval_size - 1;
  interval_fix_end = interval_end - overlap_size;

  if (ord(window) = card(window),
    interval_end     = card(tt);
    interval_fix_end = interval_end;
  );

  window_t(window,tt)     = interval_start <= ord(tt) and interval_end >= ord(tt);
  window_t_fix(window_t(window,tt)) = interval_fix_end >= ord(tt);
  t(tt)     = window_t(window,tt);
  t_fix(tt) = window_t_fix(window,tt);
  display t, t_fix;
* ***end of subset generation for rolling-horizon***


* adjusting capacities to rolling-horizon intervals
  total_plant_cap(rp)          = (%TO% - %FROM%) * yearly_plant_cap(rp) / %NBINTERVALSRH% ; // * card(t_fix)/card(t);
  total_emission_cap(e)       = (%TO% - %FROM%) * yearly_emission_cap(e) / %NBINTERVALSRH% ; //* card(t_fix)/card(t);

  Solve simple min OBJ use lp;

* Fixing time dependent variables after solve statement of timeModelWindow n
$IFTHENI.method %METHOD%==lagrange_relaxation
  FLOW.fx(t_fix,netx)          = FLOW.l(t_fix,netx)         ;
$ELSE.method
  FLOW.fx(t_fix,net)           = FLOW.l(t_fix,net)          ;
$ENDIF.method
  POWER.fx(t_fix,rp)           = POWER.l(t_fix,rp)          ;
  SLACK.fx(t_fix,r)            = SLACK.l(t_fix,r)           ;
  STORAGE_LEVEL.fx(t_fix,rs)   = STORAGE_LEVEL.l(t_fix,rs)  ;
  STORAGE_INFLOW.fx(t_fix,rs)  = STORAGE_INFLOW.l(t_fix,rs) ;
  STORAGE_OUTFLOW.fx(t_fix,rs) = STORAGE_OUTFLOW.l(t_fix,rs);

* Setting lower bounds to added capacities after solve statment of timeModelWindow n
  PLANT_ADD_CAP.lo(rp)         = PLANT_ADD_CAP.l(rp)        ;
  STORAGE_ADD_CAP.lo(rs)       = STORAGE_ADD_CAP.l(rs)      ;
$IFTHENI.method %METHOD%==lagrange_relaxation
  LINK_ADD_CAP.lo(netx)        = LINK_ADD_CAP.l(netx)       ;
$ELSE.method
  LINK_ADD_CAP.lo(net)         = LINK_ADD_CAP.l(net)        ;
$ENDIF.method
);

total_plant_cap(rp)        = (%TO% - %FROM%) * yearly_plant_cap(rp);
total_emission_cap(e)     = (%TO% - %FROM%) * yearly_emission_cap(e);
t(tt) = yes;

* If RH is called to povide good lambda for lagrange
* a) we save marginals used for lambda later BEFORE we solve the fixed LP which 'destroys' the marginals
$IFTHENI.method %METHOD%==lagrange_relaxation
  parameter rh_lambda1(e), rh_lambda2(tt,rr1,rr2), rh_lambda3(rr1,rr2), rh_obj;
  rh_lambda1(e)          = -eq_emission_cap.m(e);
  rh_lambda2(tt,rr1,rr2) = -eq_same_flow.m(tt,rr1,rr2);
  rh_lambda3(rr1,rr2)    = -eq_same_add_cap.m(rr1,rr2);
$ENDIF.method

* fixed LP is always solved
  Solve simple min OBJ use lp;

* b) we save the objective value of the fixed LP as upper bound
* c) we free fixed variables and reset bounds
$IFTHENI.method %METHOD%==lagrange_relaxation
  rh_obj = OBJ.l;
  option clear=FLOW, clear=POWER, clear=SLACK, clear=STORAGE_LEVEL, clear=STORAGE_INFLOW;
  option clear=STORAGE_OUTFLOW, clear=PLANT_ADD_CAP, clear=STORAGE_ADD_CAP, clear=LINK_ADD_CAP;
  LINK_ADD_CAP.up(netx(rr,net)) = link_max_add_cap(net)  ;
  PLANT_ADD_CAP.up(rp)          = plant_max_add_cap(rp)  ;
  STORAGE_ADD_CAP.up(rs)        = storage_max_add_cap(rs);
  STORAGE_INFLOW.up(tt,rs)      = storage_max_in(rs)     ;
  STORAGE_OUTFLOW.up(tt,rs)     = storage_max_out(rs)    ;
$ENDIF.method
