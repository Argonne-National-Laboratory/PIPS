* *** START *** Model for Benders Decomposition for Stochastic Optimization

* Based on asyncronous trust region method for stochastic programming (Linderoth 2003)

$if  not set GREENFIELD           $ set GREENFIELD        1
$if  not set USECPLEXOPT          $ set USECPLEXOPT       1
$if  not set LOADFROMXLS          $ set LOADFROMXLS       0
$if  not set USESTARINGPOINT      $ set USESTARINGPOINT   1
$if  not set USESCENARIOCUTS      $ set USESCENARIOCUTS   1
$if  not set DYNAMICTRUSTREG      $ set DYNAMICTRUSTREG   1
$if  not set DEBUGMODE            $ set DEBUGMODE         0

$onempty

Scalar loglevel;
loglevel = 2;
* Descriptions loglevel
* 0 - No display statements except starting point and final solution
* 1 - Display statements about performance of alogrithm (convergence gap)
* 2 - Display statements about points, incumbents and changes in trust region
* 3 - Display statements about submitted and running grid handles
* 4 - Display all statements

* Introduce stochastic
Set
    scen    / scen001*scen100 /
    cl      / cl001*cl001 /
    clscen(cl,scen)
    runscen(scen)
    ;

* Clustering (Clusters used in submission process as well as for cluster optimalitycuts if usescenariocuts = 0)
    loop(cl,
        loop(scen,
            if (ord(scen) > (ord(cl) - 1) * ( card(scen) /  card(cl) ) AND ord(scen) <= ord(cl) * ( card(scen) /  card(cl) ),
                clscen(cl,scen) = yes;
            );
        );
    );

* Stochastic parameters - random multiplication values for fuel prices (select one by commenting out the others and resizing the scen SET)
Parameters
* scen_cost_power_gen(scen,type,r)
* prob(scen)            / scen001 1.0 /
* cost_gas_multi(scen)  / scen001 1.0 /
* cost_coal_multi(scen) / scen001 1.0 /

 prob(scen)                       / scen001 0.2, scen002 0.4, scen003 0.2, scen004 0.15, scen005 0.05 /
 scen_type_multX(type,scen)       / gas. ( scen001 0.8, scen002 1.0, scen003 1.2, scen004 1.50, scen005 2.00 )
                                   coal. ( scen001 1.0, scen002 1.0, scen003 1.1, scen004 1.20, scen005 1.30 ) /
 scen_type_mult(scen,type)
 ;
 option scen_type_mult<scen_type_multX;

* cost_gas_multi(scen)  / scen001 0.8, scen002 1.0, scen003 1.2, scen004 1.50, scen005 2.00 /
* cost_coal_multi(scen) / scen001 1.0, scen002 1.0, scen003 1.1, scen004 1.20, scen005 1.30 /

* prob(scen)             / scen001 0.015, scen002 0.030, scen003 0.060, scen004 0.030, scen005 0.015, scen006 0.035, scen007 0.070, scen008 0.140, scen009 0.070, scen010 0.035, scen011 0.025, scen012 0.050, scen013 0.100, scen014 0.050, scen015 0.025, scen016 0.015, scen017 0.030, scen018 0.060, scen019 0.030, scen020 0.015, scen021 0.010, scen022 0.020, scen023 0.040, scen024 0.020, scen025 0.010 /
* cost_gas_multi(scen)   / scen001 0.500, scen002 0.750, scen003 1.000, scen004 1.250, scen005 1.500, scen006 0.500, scen007 0.750, scen008 1.000, scen009 1.250, scen010 1.500, scen011 0.500, scen012 0.750, scen013 1.000, scen014 1.250, scen015 1.500, scen016 0.500, scen017 0.750, scen018 1.000, scen019 1.250, scen020 1.500, scen021 0.500, scen022 0.750, scen023 1.000, scen024 1.250, scen025 1.500 /
* cost_coal_multi(scen)  / scen001 0.500, scen002 0.500, scen003 0.500, scen004 0.500, scen005 0.500, scen006 0.750, scen007 0.750, scen008 0.750, scen009 0.750, scen010 0.750, scen011 1.000, scen012 1.000, scen013 1.000, scen014 1.000, scen015 1.000, scen016 1.250, scen017 1.250, scen018 1.250, scen019 1.250, scen020 1.250, scen021 1.500, scen022 1.500, scen023 1.500, scen024 1.500, scen025 1.500 /

*prob(scen)              / scen001       0.0100  , scen002       0.0100  , scen003       0.0100  , scen004       0.0100  , scen005       0.0100  , scen006       0.0100  , scen007       0.0100  , scen008       0.0100  , scen009       0.0100  , scen010       0.0100  , scen011       0.0100  , scen012       0.0100  , scen013       0.0100  , scen014       0.0100  , scen015       0.0100  , scen016       0.0100  , scen017       0.0100  , scen018       0.0100  , scen019       0.0100  , scen020       0.0100  , scen021       0.0100  , scen022       0.0100  , scen023       0.0100  , scen024       0.0100  , scen025       0.0100  , scen026       0.0100  , scen027       0.0100  , scen028       0.0100  , scen029       0.0100  , scen030       0.0100  , scen031       0.0100  , scen032       0.0100  , scen033       0.0100  , scen034       0.0100  , scen035       0.0100  , scen036       0.0100  , scen037       0.0100  , scen038       0.0100  , scen039       0.0100  , scen040       0.0100  , scen041       0.0100  , scen042       0.0100  , scen043       0.0100  , scen044       0.0100  , scen045       0.0100  , scen046       0.0100  , scen047       0.0100  , scen048       0.0100  , scen049       0.0100  , scen050       0.0100  , scen051       0.0100  , scen052       0.0100  , scen053       0.0100  , scen054       0.0100  , scen055       0.0100  , scen056       0.0100  , scen057       0.0100  , scen058       0.0100  , scen059       0.0100  , scen060       0.0100  , scen061       0.0100  , scen062       0.0100  , scen063       0.0100  , scen064       0.0100  , scen065       0.0100  , scen066       0.0100  , scen067       0.0100  , scen068       0.0100  , scen069       0.0100  , scen070       0.0100  , scen071       0.0100  , scen072       0.0100  , scen073       0.0100  , scen074       0.0100  , scen075       0.0100  , scen076       0.0100  , scen077       0.0100  , scen078       0.0100  , scen079       0.0100  , scen080       0.0100  , scen081       0.0100  , scen082       0.0100  , scen083       0.0100  , scen084       0.0100  , scen085       0.0100  , scen086       0.0100  , scen087       0.0100  , scen088       0.0100  , scen089       0.0100  , scen090       0.0100  , scen091       0.0100  , scen092       0.0100  , scen093       0.0100  , scen094       0.0100  , scen095       0.0100  , scen096       0.0100  , scen097       0.0100  , scen098       0.0100  , scen099       0.0100  , scen100       0.0100  /
*cost_gas_multi(scen)    / scen001       0.6140  , scen002       1.8220  , scen003       1.8530  , scen004       0.7200  , scen005       1.0180  , scen006       0.7200  , scen007       1.2050  , scen008       1.6220  , scen009       1.3210  , scen010       2.1890  , scen011       1.5550  , scen012       1.0540  , scen013       2.3530  , scen014       1.5020  , scen015       2.5800  , scen016       1.8000  , scen017       0.8110  , scen018       1.8790  , scen019       2.5750  , scen020       2.2880  , scen021       1.9060  , scen022       2.6080  , scen023       2.3300  , scen024       2.9820  , scen025       0.7620  , scen026       2.8520  , scen027       2.2660  , scen028       1.0120  , scen029       2.1740  , scen030       1.6770  , scen031       1.4550  , scen032       1.0420  , scen033       1.2000  , scen034       2.3000  , scen035       1.4000  , scen036       0.9870  , scen037       2.8050  , scen038       2.9550  , scen039       1.5790  , scen040       0.6270  , scen041       1.5070  , scen042       1.8130  , scen043       2.6840  , scen044       1.0110  , scen045       1.6950  , scen046       0.9740  , scen047       1.2280  , scen048       2.6530  , scen049       1.2010  , scen050       2.5680  , scen051       2.1520  , scen052       0.7310  , scen053       2.2760  , scen054       2.2260  , scen055       2.7570  , scen056       1.8530  , scen057       1.1140  , scen058       1.0240  , scen059       1.9700  , scen060       1.9570  , scen061       1.8150  , scen062       2.2790  , scen063       2.8660  , scen064       1.4630  , scen065       1.1160  , scen066       2.1920  , scen067       2.9350  , scen068       1.4110  , scen069       0.8070  , scen070       2.8570  , scen071       1.3400  , scen072       0.9180  , scen073       1.7590  , scen074       0.7760  , scen075       0.6110  , scen076       2.2480  , scen077       1.7850  , scen078       2.5860  , scen079       1.0770  , scen080       2.7750  , scen081       1.4270  , scen082       2.3250  , scen083       2.5790  , scen084       1.2390  , scen085       2.8460  , scen086       2.5940  , scen087       1.4980  , scen088       1.8790  , scen089       0.7780  , scen090       1.4450  , scen091       1.7530  , scen092       2.7590  , scen093       1.2090  , scen094       2.3950  , scen095       1.4190  , scen096       2.1640  , scen097       0.8880  , scen098       2.6770  , scen099       1.5540  , scen100       1.6770  /
*cost_coal_multi(scen)   / scen001       0.9910  , scen002       1.0140  , scen003       1.6100  , scen004       2.6430  , scen005       2.0500  , scen006       0.9370  , scen007       1.3000  , scen008       1.1010  , scen009       2.3920  , scen010       0.6100  , scen011       1.5560  , scen012       1.0550  , scen013       2.7040  , scen014       2.6480  , scen015       2.9180  , scen016       1.8820  , scen017       1.0830  , scen018       0.7930  , scen019       0.8580  , scen020       0.9590  , scen021       2.3170  , scen022       1.1280  , scen023       1.0060  , scen024       2.0230  , scen025       1.8500  , scen026       2.7760  , scen027       1.9680  , scen028       0.7380  , scen029       2.2510  , scen030       1.9380  , scen031       2.6860  , scen032       2.1360  , scen033       0.9030  , scen034       1.1620  , scen035       2.1920  , scen036       1.7490  , scen037       2.2830  , scen038       0.7260  , scen039       1.5270  , scen040       2.7960  , scen041       0.7240  , scen042       1.8650  , scen043       1.8880  , scen044       2.3560  , scen045       2.4810  , scen046       2.9040  , scen047       2.9110  , scen048       2.8260  , scen049       2.4360  , scen050       0.7260  , scen051       2.6470  , scen052       1.1660  , scen053       0.9660  , scen054       2.8260  , scen055       1.9950  , scen056       1.8280  , scen057       2.8190  , scen058       2.2120  , scen059       0.8740  , scen060       2.6890  , scen061       1.9630  , scen062       1.0290  , scen063       1.2140  , scen064       2.1610  , scen065       1.3370  , scen066       1.4010  , scen067       2.7080  , scen068       2.7110  , scen069       0.9060  , scen070       1.2430  , scen071       0.9850  , scen072       0.7270  , scen073       1.5800  , scen074       2.1820  , scen075       2.0940  , scen076       2.0120  , scen077       1.3560  , scen078       1.8090  , scen079       2.2950  , scen080       2.4160  , scen081       2.3560  , scen082       2.0130  , scen083       2.3940  , scen084       2.8850  , scen085       2.5630  , scen086       0.6910  , scen087       1.5300  , scen088       1.9130  , scen089       2.0650  , scen090       0.9990  , scen091       2.2190  , scen092       1.1800  , scen093       0.8530  , scen094       2.3620  , scen095       1.3480  , scen096       1.4870  , scen097       1.3950  , scen098       1.5310  , scen099       0.9520  , scen100       1.7490  /
;

abort$(abs(1 - sum(scen, prob(scen))) > 1E-8) 'probabilities do not add up';


Set k                           'maximum iterations'                    / k1*k300 /;
    Alias(k,kk,k1);
Set
    k_act(k)                    'active iterations'                     / /
    k_next(k)                   'next point suggested'                  / /
    k_curinc(k)                 'current incumbent'                     / /

    k_incmap(k,k)               'incumbent for current iteration'       / /
    cur_incmap(k)               'current slice of inc_iter mapping'     / /

    kcl_done(k,cl)              'iteration-cluster finished'            / /
    cur_kcl(cl)                 'current slice of kcl mapping'          / /

    kcl_cut(k,cl)               'active cluster cuts'                   / /
    kscen_cut(k,scen)           'active scenario cuts'                  / /
    ;
Parameters
    h_k(k)                      'handles masterprob'                    / /
    h_kcl(k,cl)                 'handles subprob'                       / /

    bconst(k,scen)              'benders const'                         / /

    bcoef_plant(k,scen,rr,p)    'benders coeff'                         / /
    bcoef_stor(k,scen,rr,s)     'benders coeff'                         / /
    bcoef_link(k,scen,rr1,rr2)  'benders coeff'                         / /
    ;
$offempty

Variables
    zmaster
    zsub
    ;
Positive Variables
    theta(scen)
    ;

* Parameter for trust region algorithm
Parameters
    done            'binary indicating model is converged               '
    sub_submitted   'binary indicating sub has been submitted           '
    sub_collected   'binary indicating sub has been collected           '
    k_act_max       'maximum number of points evaluated simultaneously  '
    trr(k)          'trust region radius                                '
    trr_cur         'current trust region radius                        '
    trr_max         'maximum trust region radius                        '
    trr_min         'minimum trust region radius                        '
    trr_counter     'couter for trust region decrease                   '
    trr_rho         'evaluation for trust region                        '
    eta             'necessary improvement for new major                '
    alpha           'necessary scen fraction solved for new evaluation  '
    t_e             'time elapsed after full running                    '
    t_e_k(k)        'time elapsed at full evaluation of point           '
    ;


* Parameters for each iteration (decision variables, objective performance )
Parameters
    p_PLANT_ADD_CAP(k,rr,p),
    p_STORAGE_ADD_CAP(k,rr,s),
    p_LINK_ADD_CAP(k,rr1,rr2),
    obj_total(k),
    obj_master(k),
    obj_sub(k,scen),
    obj_model(k,scen),
    conv_gap(k)
    ;

Equations
*    Eq_objmasterrobj
    Eq_objmaster
    Eq_objsub
    Eq_objmodel
    Eq_optcut_cl(k,cl)
    Eq_optcut_scen(k,scen)
    Eq_trustregion_plant_ub(k,rr,p)
    Eq_trustregion_plant_lb(k,rr,p)
    Eq_trustregion_storage_ub(k,rr,s)
    Eq_trustregion_storage_lb(k,rr,s)
    Eq_trustregion_link_ub(k,rr1,rr2)
    Eq_trustregion_link_lb(k,rr1,rr2)
    ;

* Objective function master: Optimize generation capacities

*eq_masterrobj(r)..
*    ROBJ(r) =e= sum(rp(r,p),     PLANT_ADD_CAP(rp)   * cost_plant_add(rp))
*              + sum(rs(r,s),     STORAGE_ADD_CAP(rs) * cost_storage_add(rs));

Eq_objmaster ..
    zmaster =e= sum(rp,   PLANT_ADD_CAP(rp)   * cost_plant_add(rp))
              + sum(rs,   STORAGE_ADD_CAP(rs) * cost_storage_add(rs))
              + sum(net,  LINK_ADD_CAP(net)   * cost_link_add(net))
              + sum(scen, prob(scen) * theta(scen));

*    zmaster =e= sum((type,r),  PLANT_ADD_CAP_mod(type,r)   * cost_plant_add_mod(type,r))
*              + sum(r,  STORAGE_ADD_CAP_mod(r) * cost_storage_add_mod(r))
*              + sum(rr, LINK_ADD_CAP_mod(rr)   * cost_link_add_mod(rr))
*              + sum(scen, prob(scen) * theta(scen));

* Trust region constraints
Eq_trustregion_plant_ub(k_curinc,rp) ..
    PLANT_ADD_CAP(rp) =l= p_PLANT_ADD_CAP(k_curinc,rp) + trr_cur;
Eq_trustregion_plant_lb(k_curinc,rp) ..
    PLANT_ADD_CAP(rp) =g= p_PLANT_ADD_CAP(k_curinc,rp) - trr_cur;
Eq_trustregion_storage_ub(k_curinc,rs) ..
    STORAGE_ADD_CAP(rs) =l= p_STORAGE_ADD_CAP(k_curinc,rs) + trr_cur;
Eq_trustregion_storage_lb(k_curinc,rs) ..
    STORAGE_ADD_CAP(rs) =g= p_STORAGE_ADD_CAP(k_curinc,rs) - trr_cur;
Eq_trustregion_link_ub(k_curinc,net) ..
    LINK_ADD_CAP(net) =l= p_LINK_ADD_CAP(k_curinc,net) + trr_cur;
Eq_trustregion_link_lb(k_curinc,net) ..
    LINK_ADD_CAP(net) =g= p_LINK_ADD_CAP(k_curinc,net) - trr_cur;

* Objective function sub: Optimize dispatch
Eq_objsub ..

    zsub =e= sum((t,ptype(rp(r,p),type)), POWER(t,rp) * cost_power_generation(rp) * type_mult(type))
           + sum((t,r),                   SLACK(t,r)  * cost_unserved_demand(t))
           + sum((r,e),                   EMISSION_COST(r,e));

*    zsub =e= sum((tt,type,r), POWER_mod(tt,type,r)   * cost_power_generation_mod(type,r))
*           + sum((r,e), EMISSION(r,e)   * cost_emission(e))
*           + sum((tt,r), SLACK(tt,r)    * cost_unserved_demand(tt))
*                  ;

* Benders optimality cuts (Clustered multicut)
Eq_optcut_cl(k,cl)$kcl_cut(k,cl) ..
    sum(scen$clscen(cl,scen), theta(scen)) =g= sum(scen$clscen(cl,scen), bconst(k,scen) + sum( rp, bcoef_plant(k,scen,rp) * PLANT_ADD_CAP(rp) )
                                                                                        + sum( rs, bcoef_stor(k,scen,rs) * STORAGE_ADD_CAP(rs) )
                                                                                        + sum( net, bcoef_link(k,scen,net) * LINK_ADD_CAP(net) ) );

* Benders optimality cuts (Scenario multicut)
Eq_optcut_scen(k,scen)$kscen_cut(k,scen) ..
    theta(scen) =g= bconst(k,scen) + sum( rp, bcoef_plant(k,scen,rp) * PLANT_ADD_CAP(rp) )
                                   + sum( rs, bcoef_stor(k,scen,rs) * STORAGE_ADD_CAP(rs) )
                                   + sum( net, bcoef_link(k,scen,net) * LINK_ADD_CAP(net) );


model master / Eq_objmaster, Eq_optcut_cl, Eq_optcut_scen, Eq_trustregion_plant_ub, Eq_trustregion_plant_lb, Eq_trustregion_storage_ub, Eq_trustregion_storage_lb, Eq_trustregion_link_ub, Eq_trustregion_link_lb /;
model sub / simple - eq_obj - eq_robj + Eq_objsub /;

* Options for solver
option limrow=0, limcol=0, solprint=silent;

* No option for asyncThreads due to compatibility issues with solvelink.AsyncThreads and GUSS
option lp = cplex;
simple.solvelink = %solvelink.ChainScript%;
master.solvelink = %solvelink.AsyncGrid%;
sub.solvelink = %solvelink.AsyncGrid%;


* Effect of holdfixed on model building time?
simple.holdfixed = 1;
sub.holdfixed = 1;
$IFTHEN.useCPLEXopt %USECPLEXOPT%==1
simple.optfile = 1;
master.optfile = 1;
sub.optfile = 1;
$ENDIF.useCPLEXopt

* *** END *** Model for Benders Decomposition


$IFTHEN.useSTARTINGPOINT %USESTARINGPOINT%==1
* Region independent calculation of starting point
LINK_ADD_CAP.fx(net) = 0;
solve simple min OBJ us lp;
LINK_ADD_CAP.lo(net) = 0;
LINK_ADD_CAP.up(net) = link_max_add_cap(net);

display "-----------------------------------------------------------------",
        "Starting point for independent regions",
        "-----------------------------------------------------------------",
        PLANT_ADD_CAP.l,STORAGE_ADD_CAP.l,LINK_ADD_CAP.l;
* Initialization of major starting point
p_PLANT_ADD_CAP('k1',rp)   = PLANT_ADD_CAP.l(rp);
p_STORAGE_ADD_CAP('k1',rs) = STORAGE_ADD_CAP.l(rs)   ;
p_LINK_ADD_CAP('k1',net)   = LINK_ADD_CAP.l(net)     ;
$ELSE.useSTARTINGPOINT

display "-----------------------------------------------------------------",
        "Using trivial starting point",
        "-----------------------------------------------------------------",
        PLANT_ADD_CAP.l,STORAGE_ADD_CAP.l,LINK_ADD_CAP.l;
* Initialization of major starting point
p_PLANT_ADD_CAP('k1',rp)   = 0;
p_STORAGE_ADD_CAP('k1',rs) = 0;
p_LINK_ADD_CAP('k1',net)   = 0;
$ENDIF.useSTARTINGPOINT


* ************************
* Parameters and configuration for GUSS usage
* ************************
Parameters
*s_cost_power_gen(scen,type,r)
*s_POWER_mod(scen,t,type,r)
*s_STORAGE_LEVEL_mod(scen,t,r)
*s_STORAGE_INFLOW_mod(scen,t,r)
*s_STORAGE_OUTFLOW_mod(scen,t,r)
*s_FLOW_mod(scen,t,r1,r2)
*s_SLACK(scen,t,r)
*s_EMISSION(scen,r,e)

s_zsub(scen)
s_eq_power_balance(scen,tt,rr)
s_eq_plant_capacity(scen,tt,rr,p)
s_eq_total_plant_capacity(scen,rr,p)
s_eq_storage_capacity(scen,tt,rr,s)
s_eq_emission_cap(scen,e)
s_eq_link_capacity(scen,tt,rr1,rr2)
s_eq_storage_bound_inflow(scen,tt,rr,s)
s_eq_storage_bound_outflow(scen,tt,rr,s)
;
*s_cost_power_gen(scen,type,r) = scen_cost_power_gen(scen,type,r);

Parameter
g_opt / SkipBaseCase 0, OptfileInit 1 /

* Set of scenarios and GUSS options
SET dict /  runscen                     .scenario   .''
            g_opt                       .opt        .''
* Updater for stochastic parameters
            type_mult                   . param     . scen_type_mult

* Gatherer for sub scenario variables (maybe can be skipped to improve GUSS speed)
*            POWER_mod                   .level      .s_POWER_mod
*            STORAGE_LEVEL_mod           .level      .s_STORAGE_LEVEL_mod
*            STORAGE_INFLOW_mod          .level      .s_STORAGE_INFLOW_mod
*            STORAGE_OUTFLOW_mod         .level      .s_STORAGE_OUTFLOW_mod
*            FLOW_mod                    .level      .s_FLOW_mod
*            SLACK                       .level      .s_SLACK
*            EMISSION                    .level      .s_EMISSION

* Gatherer for sub_obj level and marginals (maybe can be reduced to marginals of first stage vars to improve GUSS speed)
            zsub                        .level      .s_zsub
            eq_plant_capacity           .marginal   .s_eq_plant_capacity
            eq_storage_capacity         .marginal   .s_eq_storage_capacity
            eq_link_capacity            .marginal   .s_eq_link_capacity

*            eq_power_balance            .marginal   .s_eq_power_balance
*            eq_total_plant_capacity     .marginal   .s_eq_total_plant_capacity
*            eq_emission_cap             .marginal   .s_eq_emission_cap
*            eq_storage_bound_inflow     .marginal   .s_eq_storage_bound_inflow
*            eq_storage_bound_outflow    .marginal   .s_eq_storage_bound_outflow
/;
display$(loglevel >= 4) dict;


* ************************
* Configuration of Benders algorithm
* ************************
* Set finishing objective = no
done = 0;
* maximum trust region radius
trr_counter = 0;
* Set starting incumbent
k_curinc('k1') = yes;
* Set staring active iteration range
k_act('k1') = yes;


* Allowed discrepancy between Q and m in order to accept new point
eta = 0.01;
* Maximum number of points beeing evaluated simultaneously
k_act_max = 3;
* neccessary fraction solved in order to start new evaluation
alpha = 0.4;

* initial trust region radius ('k1' equals starting point)
trr(k) = 1;
trr('k1') = 0;
* limits on dynamic trust-region size
trr_max = 100;
trr_min = 0.01;


* ************************
* Benders algorithm
* ************************

* Initialize necessary parameters for syntax (values will be overwritten)
obj_sub('k1',scen) = 0;

* Loop until done or noch active iterations left
while( card(k_act) > 0 AND NOT done ,

* ************************
* Masterproblems Submit
* ************************

* Loop all currently active iterations
    loop(kk$(k_act(kk) AND NOT done),

* Get slice of incumbent map for current iteration
        cur_incmap(k1) = k_incmap(k1,kk);
* If slice of incumbent map is empty, point has not been suggested yet
        if(card(cur_incmap) = 0,

* Set current incumbent for current iteration
            k_incmap(k_curinc,kk) = yes;
            display$(loglevel >= 2) 'Current incumbent mapped',
                    k_incmap;

* Get current trust region
            trr_cur = trr(kk);

* Solve master problem to propose x(kk) consider trust region of current incumbent (k_curinc)
            solve master min zmaster us lp;

            h_k(kk) = master.handle;
            display$(loglevel >= 4) 'New master handle submitted',
                    h_k;
        );
    );

* ************************
* Masterproblems Collect
* ************************

* Loop all currently active iterations
    loop(kk$handlecollect(h_k(kk)),

        display$(loglevel >= 4) 'Master handle collected',
            h_k;
        display$handledelete(h_k(kk)) 'ERROR: Trouble deleting master handle' ;
        h_k(kk) = 0;

        obj_master(kk) = zmaster.l - sum(scen, prob(scen) * theta.l(scen));
        obj_model(kk,scen) = theta.l(scen);

* Check convergence (Check if starting point has been fully evaluated yet)
        cur_kcl(cl) = kcl_done('k1',cl);
        if(card(cur_kcl) = card(cl),
            done = sum( (scen,k1)$k_curinc(k1), prob(scen) * obj_sub(k1,scen) ) + sum( k1$k_curinc(k1), obj_master(k1) ) - (sum(scen, prob(scen) * obj_model(kk,scen)) + obj_master(kk)) <= 1e-8 * ( 1 + sum( (scen,k1)$k_curinc(k1), prob(scen) * obj_sub(k1,scen) ) + sum( k1$k_curinc(k1), obj_master(k1) ) );
        );

* Save current solution
        p_PLANT_ADD_CAP(kk,rp)    = PLANT_ADD_CAP.l(rp)  ;
        p_STORAGE_ADD_CAP(kk,rs)  = STORAGE_ADD_CAP.l(rs);
        p_LINK_ADD_CAP(kk,net)    = LINK_ADD_CAP.l(net)  ;
    );

* ************************
* Subproblems Submit
* ************************
* Loop all currently active iterations
* (one sub collected per main loop)

    sub_submitted = 0;
    loop(kk$(k_act(kk) AND NOT done AND NOT sub_submitted),

* Get slice of incumbent map for current iteration
        cur_incmap(k1) = k_incmap(k1,kk);
* Check if point has incumbent and handle for masterproblem is empty (masterproblem has been solved)
        if(card(cur_incmap) = 1 AND h_k(kk) = 0,

* Fix bounds
            PLANT_ADD_CAP.fx(rp)    = p_PLANT_ADD_CAP(kk,rp) ;
            STORAGE_ADD_CAP.fx(rs)  = p_STORAGE_ADD_CAP(kk,rs);
            LINK_ADD_CAP.fx(net)    = p_LINK_ADD_CAP(kk,net);

* Submit not solved scenario clusters
            loop(cl$( h_kcl(kk,cl) = 0 AND NOT kcl_done(kk,cl) AND NOT sub_submitted ),
                runscen(scen) = clscen(cl,scen);

* Evaluate subgradient at x(kk) for each scenario
                solve sub min zsub us lp scenario dict;

                h_kcl(kk,cl) = sub.handle;
                display$(loglevel >= 4) 'New sub handle submitted',
                        h_kcl;

                sub_submitted = 1;
                );
        );

* Reset bounds
            PLANT_ADD_CAP.lo(rp)    = 0;
            STORAGE_ADD_CAP.lo(rs)  = 0;
            LINK_ADD_CAP.lo(net)    = 0;
            PLANT_ADD_CAP.up(rp)    = plant_max_add_cap(rp);
            STORAGE_ADD_CAP.up(rs)  = storage_max_add_cap(rs);
            LINK_ADD_CAP.up(net)    = link_max_add_cap(net);
    );

* ************************
* Subproblems Collect
* ************************
* Loop all currently active iterations, sub_collected forces to skip current collection loop
* (one sub collected per main loop - GUSS scattering is computational bottleneck and prevents new points to be evaluated before the old one)

    sub_collected = 0;
    loop(kk$(k_act(kk) AND NOT done AND NOT sub_collected),

        loop(cl$(handlecollect(h_kcl(kk,cl)) AND NOT sub_collected),


            display$(loglevel >= 4) 'Sub handle collected',
                                    h_kcl;
            display$handledelete(h_kcl(kk,cl)) 'ERROR: Trouble deleting sub handle' ;
                h_kcl(kk,cl) = 0;
            sub_collected = 1;

            obj_sub(kk,scen) = s_zsub(scen);

            bcoef_plant(kk,scen,rp)  = sum(t, s_eq_plant_capacity(scen,t,rp) );
            bcoef_stor(kk,scen,rs)   = sum(t, s_eq_storage_capacity(scen,t,rs) );
            bcoef_link(kk,scen,net)  = sum(t, s_eq_link_capacity(scen,t,net)  );

            bconst(kk,scen) =  obj_sub(kk,scen)
                             - sum(rp,  bcoef_plant(kk,scen,rp) * p_PLANT_ADD_CAP(kk,rp))
                             - sum(rs,  bcoef_stor(kk,scen,rs)  * p_STORAGE_ADD_CAP(kk,rs))
                             - sum(net, bcoef_link(kk,scen,net) * p_LINK_ADD_CAP(kk,net));

*            bconst(kk,scen) = sum( (t,r), s_eq_power_balance(scen,t,r) * ( renewables(t,r) - demand(t,r) ) )
*                            + sum( (rp,t), s_eq_plant_capacity(scen,t,rp) * plant_cap(rp) )
*                            + sum( rp, s_eq_total_plant_capacity(scen,rp) * total_plant_cap(rp) )
*                            + sum( (rs,t), s_eq_storage_capacity(scen,t,rs)  * storage_cap_mod(rs) )
*                            + sum( e, s_eq_emission_cap(scen,e) * total_emission_cap(e) )
*                            + sum( (net,t), s_eq_link_capacity(scen,t,net) * link_cap(net) )
*                            + sum( (rs,t), s_eq_storage_bound_inflow(scen,t,rs) * storage_max_in(rs) )
*                            + sum( (rs,t), s_eq_storage_bound_outflow(scen,t,rs) * storage_max_out(rs) )
*                            ;


* Set kcl_done = yes for all evaluated clusters in iteration kk
            kcl_done(kk,cl) = yes;

$IFTHEN.useScenarioCuts %USESCENARIOCUTS%==1
            kscen_cut(kk,scen) = yes;
$ELSE.useScenarioCuts
            kcl_cut(kk,cl) = yes;
$ENDIF.useScenarioCuts
        );
    );


* ************************
* Evaluation of point
* ************************
* Loop all currently active iterations
    loop(kk$(k_act(kk) AND NOT done),

* ************************
* Partial evaluation
* ************************
* Get slice of kcl_done mapping
        cur_kcl(cl) = kcl_done(kk,cl);
* If sufficient clusters have been evaluated activate the following iteratation
        if ( card(cur_kcl) >= alpha * card(cl) AND card(k_act) < k_act_max AND NOT k_next(kk),
            display$(loglevel >= 2) 'Point partially evaluated, starting next point',
                    cur_kcl;
            k_act(kk+1) = yes;
            k_next(kk) = yes;
            display$(loglevel >= 3) 'New set of active iterations',
                    k_act;
        );


* ************************
* Full evaluation
* ************************
* Get slice of kcl_done mapping
        cur_kcl(cl) = kcl_done(kk,cl);
* If all clusters have been evaluated...
        if (card(cur_kcl) = card(cl),
            display$(loglevel >= 2) 'Point fully evaluated',
                    cur_kcl;
            obj_total(kk) = obj_master(kk) + sum(scen, prob(scen) * obj_sub(kk,scen));
            display$(loglevel >= 4) obj_total, obj_master, obj_model, obj_sub;

* ...else add optimality cut to model and possibly reduce trust region radius
            if( sum(scen, prob(scen) * obj_sub(kk,scen)) + obj_master(kk)  >  sum( (scen,k1)$k_curinc(k1), prob(scen) * obj_sub(k1,scen) ) + sum( k1$k_curinc(k1), obj_master(k1) ) - eta * ( (sum( (scen,k1)$k_curinc(k1), prob(scen) * obj_sub(k1,scen) ) + sum( k1$k_curinc(k1), obj_master(k1) )) - (sum(scen, prob(scen) * obj_model(kk,scen)) + obj_master(kk)) ) ,

* Decrease trust region for large difference between model and sub objective or consecutive errors
$IFTHEN.useDynamicTrustRegion %USEDYNAMICTRUSTREG%==1
                trr_rho = min(1,trr(k)) * ( (sum(scen, prob(scen)*obj_sub(kk,scen))+obj_master(kk)) - (sum((scen,k1)$k_curinc(k1),prob(scen)*obj_sub(k1,scen))+sum(k1$k_curinc(k1),obj_master(k1))) ) / ( (sum((scen,k1)$k_curinc(k1),prob(scen)*obj_sub(k1,scen))+sum(k1$k_curinc(k1),obj_master(k1))) - (sum(scen, prob(scen)*obj_model(kk,scen))+obj_master(kk)) );

                if(trr_rho > 0,
                    trr_counter = trr_counter + 1;
                );
                if(trr_rho > 3 OR ( trr_counter >= 3 AND (trr_rho > 1 AND trr_rho <= 3) ),

                    trr(k)$( NOT k_next AND NOT k_act ) = max(trr_min, 1/min(trr_rho, 4) * trr(k));

                    display$(loglevel >= 2) 'Trust region updated',
                            trr;

                    trr_counter = 0;
                );

$ENDIF.useDynamicTrustRegion

*            );
* ...check if solution is acceped as new major and possibly increase trust region radius
* global incumbent
*            if( sum(scen, prob(scen) * obj_sub(kk,scen)) + obj_master(kk)  <=  sum( (scen,k1)$k_curinc(k1), prob(scen) * obj_sub(k1,scen) ) + sum( k1$k_curinc(k1), obj_master(k1) ) - eta * ( (sum( (scen,k1)$k_curinc(k1), prob(scen) * obj_sub(k1,scen) ) + sum( k1$k_curinc(k1), obj_master(k1) )) - (sum(scen, prob(scen) * obj_model(kk,scen)) + obj_master(kk)) ) ,
            else
                k_curinc(k) = no;
                k_curinc(kk) = yes;

                display$(loglevel >= 2) 'New incumbent set',
                        k_curinc;

* Increase trust region if model approximation and actual sub objective are close
$IFTHEN.useDynamicTrustRegion %USEDYNAMICTRUSTREG%==1
* CHECK: k_curinc already has kk as current incumbent
                if( sum(scen, prob(scen) * obj_sub(kk,scen)) + obj_master(kk)  <=  sum( (scen,k1)$k_curinc(k1), prob(scen) * obj_sub(k1,scen) ) + sum( k1$k_curinc(k1), obj_master(k1) ) - 0.5 * ( (sum( (scen,k1)$k_curinc(k1), prob(scen) * obj_sub(k1,scen) ) + sum( k1$k_curinc(k1), obj_master(k1) )) - (sum(scen, prob(scen) * obj_model(kk,scen)) + obj_master(kk)) ) ,
                    trr(k)$( NOT k_next AND NOT k_act ) = min( 2 * trr(kk),trr_max );

                    display$(loglevel >= 2) 'Trust region updated',
                            trr;

                );
$ENDIF.useDynamicTrustRegion

            );

            k_act(kk) = no;

            conv_gap(k1)$( k_next(k1) AND NOT k_act(k1) AND sum(scen,obj_model(k1,scen)) > 0 ) = sum(scen, prob(scen) * (obj_sub(k1,scen) - obj_model(k1,scen)) ) / ( sum(scen, prob(scen) * obj_model(k1,scen)) + obj_master(k1) );
            display$(loglevel >= 1) 'Current convergence gap:',
                    conv_gap;
            display$(loglevel >= 3) 'New set of active iterations',
                    k_act;

* Log time of point evaluated
        t_e_k(kk) = timeelapsed;
        );
    );

* Waiting phase
    display$(loglevel >= 3) 'Handles currently running:',
            h_k, h_kcl;

    display$sleep(1) 'Waiting';
);


* ************************
* End of Benders iteration loop
* ************************

if(card(k_act) = 0 AND NOT done,
display "-----------------------------------------------------------------",
        "Limit for iterations reached",
        "-----------------------------------------------------------------";
);
if(card(k_act) > 0 AND done,
display "-----------------------------------------------------------------",
        "Optimal solution found",
        "-----------------------------------------------------------------";
);

* Get final solution
Parameter
final_PLANT_ADD_CAP(rr,p)
final_STORAGE_ADD_CAP(rr,s)
final_LINK_ADD_CAP(rr1,rr2)
final_obj_total
final_obj_master
final_obj_sub;

final_PLANT_ADD_CAP(rp)   = sum(k$k_curinc(k), p_PLANT_ADD_CAP(k,rp)                  );
final_STORAGE_ADD_CAP(rs) = sum(k$k_curinc(k), p_STORAGE_ADD_CAP(k,rs)                );
final_LINK_ADD_CAP(net)   = sum(k$k_curinc(k), p_LINK_ADD_CAP(k,net)                  );
final_obj_total           = sum(k$k_curinc(k), obj_total(k)                           );
final_obj_master          = sum(k$k_curinc(k), obj_master(k)                          );
final_obj_sub             = sum(k$k_curinc(k), sum(scen,prob(scen) * obj_sub(k,scen)) );

display "-----------------------------------------------------------------",
        "Trust region solution and cost",
        final_PLANT_ADD_CAP,final_STORAGE_ADD_CAP,final_LINK_ADD_CAP,final_obj_total,final_obj_master,final_obj_sub,
        "-----------------------------------------------------------------";
t_e = timeelapsed;
display "-----------------------------------------------------------------",
        " Time elapsed: ", t_e,
        "-----------------------------------------------------------------";

Execute_Unload 'atr_cl_simple_guss_results_%NBREGIONS%.gdx',done,t_e,t_e_k,trr,k_act,clscen,kcl_done,kcl_cut,kscen_cut,k_incmap,obj_total,obj_master,obj_model,obj_sub,bconst,bcoef_plant,bcoef_stor,bcoef_link,p_PLANT_ADD_CAP,p_STORAGE_ADD_CAP,p_LINK_ADD_CAP,final_PLANT_ADD_CAP,final_STORAGE_ADD_CAP,final_LINK_ADD_CAP,final_obj_total,final_obj_master,final_obj_sub;
