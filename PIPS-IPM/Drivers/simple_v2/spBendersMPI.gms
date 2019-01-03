$if  not set USECPLEXOPT          $ set USECPLEXOPT       1
$if  not set USESTARTINGPOINT     $ set USESTARTINGPOINT  1
$if  not set USETRUSTREG          $ set USETRUSTREG       1
$if  not set DYNAMICTRUSTREG      $ set DYNAMICTRUSTREG   1
$if  not set NBSCEN               $ set NBSCEN            5
$if  not set NBCLUSTER            $ set NBCLUSTER         %NBSCEN%
$if  not set BENDERSITERLIM       $ set BENDERSITERLIM    300
$if  not set CONVERGENCECRIT      $ set CONVERGENCECRIT   1E-5

$eval CLHELPER ceil(%NBSCEN%/%NBCLUSTER%)

$eval NBPROC %NBSCEN%+1
$if not %NBPROC%==%sysenv.PMI_SIZE% $abort MPI size needs to be %NBPROC%


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
    scen    / scen1*scen%NBSCEN% /
    activeScen(scen)
    cl      / cl1*cl%NBCLUSTER% /
    clscen(cl,scen);
 
Set cl_helper / clh1*clh%CLHELPER% /;
Set clscen_map(cl, cl_helper, scen)  / #cl . #cl_helper : #scen /;    
Option clscen < clscen_map;    

Parameter type_mult_scen(scen,type)

$ifthen.nbscen %NBSCEN%==5
   parameter
      prob(scen)                 / scen1 0.2, scen2 0.4, scen3 0.2, scen4 0.15, scen5 0.05 /
      scen_type_multX(type,scen) / gas. ( scen1 0.8, scen2 1.0, scen3 1.2, scen4 1.50, scen5 2.00 )
                                   coal.( scen1 1.0, scen2 1.0, scen3 1.1, scen4 1.20, scen5 1.30 )
                                   wind.#scen 1, photovoltaic.#scen 1 /
   option type_mult_scen<scen_type_multX;

$elseif.nbscen %NBSCEN%==100
   parameter
      prob(scen)
      scen_type_multX(type,scen) /
      gas. (
      scen1         0.6140  , scen2        1.8220  , scen3        1.8530  ,
      scen4         0.7200  , scen5        1.0180  , scen6        0.7200  ,
      scen7         1.2050  , scen8        1.6220  , scen9        1.3210  ,
      scen10        2.1890  , scen11       1.5550  , scen12       1.0540  ,
      scen13        2.3530  , scen14       1.5020  , scen15       2.5800  ,
      scen16        1.8000  , scen17       0.8110  , scen18       1.8790  ,
      scen19        2.5750  , scen20       2.2880  , scen21       1.9060  ,
      scen22        2.6080  , scen23       2.3300  , scen24       2.9820  ,
      scen25        0.7620  , scen26       2.8520  , scen27       2.2660  ,
      scen28        1.0120  , scen29       2.1740  , scen30       1.6770  ,
      scen31        1.4550  , scen32       1.0420  , scen33       1.2000  ,
      scen34        2.3000  , scen35       1.4000  , scen36       0.9870  ,
      scen37        2.8050  , scen38       2.9550  , scen39       1.5790  ,
      scen40        0.6270  , scen41       1.5070  , scen42       1.8130  ,
      scen43        2.6840  , scen44       1.0110  , scen45       1.6950  ,
      scen46        0.9740  , scen47       1.2280  , scen48       2.6530  ,
      scen49        1.2010  , scen50       2.5680  , scen51       2.1520  ,
      scen52        0.7310  , scen53       2.2760  , scen54       2.2260  ,
      scen55        2.7570  , scen56       1.8530  , scen57       1.1140  ,
      scen58        1.0240  , scen59       1.9700  , scen60       1.9570  ,
      scen61        1.8150  , scen62       2.2790  , scen63       2.8660  ,
      scen64        1.4630  , scen65       1.1160  , scen66       2.1920  ,
      scen67        2.9350  , scen68       1.4110  , scen69       0.8070  ,
      scen70        2.8570  , scen71       1.3400  , scen72       0.9180  ,
      scen73        1.7590  , scen74       0.7760  , scen75       0.6110  ,
      scen76        2.2480  , scen77       1.7850  , scen78       2.5860  ,
      scen79        1.0770  , scen80       2.7750  , scen81       1.4270  ,
      scen82        2.3250  , scen83       2.5790  , scen84       1.2390  ,
      scen85        2.8460  , scen86       2.5940  , scen87       1.4980  ,
      scen88        1.8790  , scen89       0.7780  , scen90       1.4450  ,
      scen91        1.7530  , scen92       2.7590  , scen93       1.2090  ,
      scen94        2.3950  , scen95       1.4190  , scen96       2.1640  ,
      scen97        0.8880  , scen98       2.6770  , scen99       1.5540  ,
      scen100       1.6770 )

      coal.(
      scen1         0.9910  , scen2        1.0140  , scen3        1.6100  ,
      scen4         2.6430  , scen5        2.0500  , scen6        0.9370  ,
      scen7         1.3000  , scen8        1.1010  , scen9        2.3920  ,
      scen10        0.6100  , scen11       1.5560  , scen12       1.0550  ,
      scen13        2.7040  , scen14       2.6480  , scen15       2.9180  ,
      scen16        1.8820  , scen17       1.0830  , scen18       0.7930  ,
      scen19        0.8580  , scen20       0.9590  , scen21       2.3170  ,
      scen22        1.1280  , scen23       1.0060  , scen24       2.0230  ,
      scen25        1.8500  , scen26       2.7760  , scen27       1.9680  ,
      scen28        0.7380  , scen29       2.2510  , scen30       1.9380  ,
      scen31        2.6860  , scen32       2.1360  , scen33       0.9030  ,
      scen34        1.1620  , scen35       2.1920  , scen36       1.7490  ,
      scen37        2.2830  , scen38       0.7260  , scen39       1.5270  ,
      scen40        2.7960  , scen41       0.7240  , scen42       1.8650  ,
      scen43        1.8880  , scen44       2.3560  , scen45       2.4810  ,
      scen46        2.9040  , scen47       2.9110  , scen48       2.8260  ,
      scen49        2.4360  , scen50       0.7260  , scen51       2.6470  ,
      scen52        1.1660  , scen53       0.9660  , scen54       2.8260  ,
      scen55        1.9950  , scen56       1.8280  , scen57       2.8190  ,
      scen58        2.2120  , scen59       0.8740  , scen60       2.6890  ,
      scen61        1.9630  , scen62       1.0290  , scen63       1.2140  ,
      scen64        2.1610  , scen65       1.3370  , scen66       1.4010  ,
      scen67        2.7080  , scen68       2.7110  , scen69       0.9060  ,
      scen70        1.2430  , scen71       0.9850  , scen72       0.7270  ,
      scen73        1.5800  , scen74       2.1820  , scen75       2.0940  ,
      scen76        2.0120  , scen77       1.3560  , scen78       1.8090  ,
      scen79        2.2950  , scen80       2.4160  , scen81       2.3560  ,
      scen82        2.0130  , scen83       2.3940  , scen84       2.8850  ,
      scen85        2.5630  , scen86       0.6910  , scen87       1.5300  ,
      scen88        1.9130  , scen89       2.0650  , scen90       0.9990  ,
      scen91        2.2190  , scen92       1.1800  , scen93       0.8530  ,
      scen94        2.3620  , scen95       1.3480  , scen96       1.4870  ,
      scen97        1.3950  , scen98       1.5310  , scen99       0.9520  ,
      scen100       1.7490 ) /
      
   option type_mult_scen<scen_type_multX;
   prob(scen) = 1/card(scen);
$else.nbscen
   parameter
      prob(scen);

* Do it for all scenarios, since we want the same random numbers in every block
   type_mult_scen(scen,type)$(not re(type)) = uniform(0.5,2.5);
   prob(scen) = 1/card(scen);

$endif.nbscen

abort$(abs(1 - sum(scen, prob(scen))) > 1E-8) 'probabilities do not add up';


Set k                           'maximum iterations'                    / k1*k%BENDERSITERLIM% /
    rephdr                      'header for report parameter'           / Obj, globObj, Bound, absGap, relGap, subTime, masterTime, Time /
    ;
    Alias(k,kk,k1);

$onempty
Set
    k_curinc(k)                 'current incumbent'                     / /
    k_incmap(k,k)               'incumbent for current iteration'       / /
    kcl_cut(k,cl)               'active cluster cuts'                   / /
    ;
Parameters
    bconst(k,scen)              'benders const'                         / /
    bconst_check(k,scen)        'benders const'                         / /
    bcoef_plant(k,scen,rr,p)    'benders coeff'                         / /
    bcoef_stor(k,scen,rr,s)     'benders coeff'                         / /
    bcoef_link(k,scen,rr1,rr2)  'benders coeff'                         / /

    sub_bconst(scen)              'benders const'                         / /
    sub_bcoef_plant(scen,rr,p)    'benders coeff'                         / /
    sub_bcoef_stor(scen,rr,s)     'benders coeff'                         / /
    sub_bcoef_link(scen,rr1,rr2)  'benders coeff'                         / /

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
    done            'binary indicating model is converged               ' /0/
    trr(k)          'trust region radius                                '
    trr_cur         'current trust region radius                        '
    trr_max         'maximum trust region radius                        '
    trr_min         'minimum trust region radius                        '
    trr_counter     'couter for trust region decrease                   '
    trr_rho         'evaluation for trust region                        '
    eta             'necessary improvement for new major                '
    alpha           'necessary scen fraction solved for new evaluation  '
    t_e             'time elapsed after full run                        '
    t_e_k(k)        'time elapsed at full evaluation of point           '
    ;


* Parameters for each iteration (decision variables, objective performance )
Parameters
    p_PLANT_ADD_CAP(k,rr,p),
    p_STORAGE_ADD_CAP(k,rr,s),
    p_LINK_ADD_CAP(k,rr1,rr2),

    sub_PLANT_ADD_CAP(rr,p),
    sub_STORAGE_ADD_CAP(rr,s),
    sub_LINK_ADD_CAP(rr1,rr2),
    
    obj_total(k),
    obj_master(k),
    obj_sub(k,scen),
    sub_obj_sub(scen),
    obj_model(k,scen),
    conv_gap(k)
    
    upperbound(k)
    lowerbound(k)
    rep(*,*) 
    ;
    
Scalar
    upperbound_glob / +INF /
    lowerbound_glob / 0 /
    ;

Equations
    Eq_objmaster
    Eq_objsub
    Eq_objmodel
    Eq_optcut(k,cl)
    Eq_trustregion_plant_ub(k,rr,p)
    Eq_trustregion_plant_lb(k,rr,p)
    Eq_trustregion_storage_ub(k,rr,s)
    Eq_trustregion_storage_lb(k,rr,s)
    Eq_trustregion_link_ub(k,rr1,rr2)
    Eq_trustregion_link_lb(k,rr1,rr2)
    ;

* Objective function master: Optimize generation capacities
Eq_objmaster ..
    zmaster =e= sum(rp,   PLANT_ADD_CAP(rp)   * cost_plant_add(rp))
              + sum(rs,   STORAGE_ADD_CAP(rs) * cost_storage_add(rs))
              + sum(net,  LINK_ADD_CAP(net)   * cost_link_add(net))
              + sum(scen, prob(scen) * theta(scen));

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

* Benders optimality cuts (Clustered multicut, Scenario cuts for nbcl = nbscen)
Eq_optcut(k,cl)$kcl_cut(k,cl) ..
    sum(scen$clscen(cl,scen), theta(scen)) =g= sum(scen$clscen(cl,scen), bconst(k,scen) + sum( rp, bcoef_plant(k,scen,rp) * PLANT_ADD_CAP(rp) )
                                                                                        + sum( rs, bcoef_stor(k,scen,rs) * STORAGE_ADD_CAP(rs) )
                                                                                        + sum( net, bcoef_link(k,scen,net) * LINK_ADD_CAP(net) ) );
                                                                          
model master / 
    Eq_objmaster, 
    Eq_optcut, 
$IFTHEN.useTrustRegion %USETRUSTREG%==1
    Eq_trustregion_plant_ub, 
    Eq_trustregion_plant_lb, 
    Eq_trustregion_storage_ub, 
    Eq_trustregion_storage_lb, 
    Eq_trustregion_link_ub, 
    Eq_trustregion_link_lb
$ENDIF.useTrustRegion
    /;
    
model sub / 
    simple 
    - eq_obj 
    - eq_robj 
    + Eq_objsub 
    /;

* Options for solver
option limrow=0, limcol=0, solprint=silent, solvelink=5;
option lp = cplexd;

simple.holdfixed = 1;
sub.holdfixed = 1;
$IFTHEN.useCPLEXopt %USECPLEXOPT%==1
$onecho > cplexd.opt
names no
*threads 4
lpmethod 4
$offecho
simple.optfile = 1;
master.optfile = 1;
sub.optfile = 1;
$ENDIF.useCPLEXopt

embeddedCode Python:
  from mpi4py import *
  comm = MPI.COMM_WORLD

  def symbols2List(*symbols):
      return list(map(lambda x: list(gams.get(x)), symbols))
pauseEmbeddedCode
file putlog / '' /;

$ifthen.MPI 0==%sysenv.PMI_RANK%
$ IFTHEN.useSTARTINGPOINT %USESTARTINGPOINT%==1
*   Use mean scenario for calculation of starting point
    type_mult(type) = sum(scen, prob(scen) * type_mult_scen(scen,type));
	putclose putlog / '##### Solve mean scenario' /;
	simple.reslim = 10;
    solve simple min OBJ us lp;
        
    display "-----------------------------------------------------------------",
            "Starting point for stochastic mean",
            "-----------------------------------------------------------------",
            PLANT_ADD_CAP.l,STORAGE_ADD_CAP.l,LINK_ADD_CAP.l;
*   Initialization of major starting point
    p_PLANT_ADD_CAP('k1',rp)   = PLANT_ADD_CAP.l(rp);
    p_STORAGE_ADD_CAP('k1',rs) = STORAGE_ADD_CAP.l(rs)   ;
    p_LINK_ADD_CAP('k1',net)   = LINK_ADD_CAP.l(net)     ;
$ ELSE.useSTARTINGPOINT
    display "-----------------------------------------------------------------",
            "Using trivial starting point",
            "-----------------------------------------------------------------",
            PLANT_ADD_CAP.l,STORAGE_ADD_CAP.l,LINK_ADD_CAP.l;
*   Initialization of major starting point
    p_PLANT_ADD_CAP('k1',rp)   = 0;
    p_STORAGE_ADD_CAP('k1',rs) = 0;
    p_LINK_ADD_CAP('k1',net)   = 0;
$ENDIF.useSTARTINGPOINT

* Configuration of Benders algorithm

* Set finishing objective = no
  done = 0;

* Set starting incumbent
  k_curinc('k1') = yes;

* Allowed discrepancy between Q and m in order to accept new point (0.01)
  eta = 0;

* Initial trust region radius ('k1' equals starting point)
  trr(k) = 1;
  trr('k1') = 0;

* Trust region counter
  trr_counter = 0;

* Limits on dynamic trust-region size
  trr_max = 10;
  trr_min = 0.1;
file fff;
* Loop all currently active iterations
loop(kk$( NOT done ),
    rep(kk,'masterTime') = timeelapsed;
* Set current incumbent for current iteration
    k_incmap(k_curinc,kk) = yes;

* Get current trust region
    trr_cur = trr(kk);
    put_utility 'gdxout' / 'mpiPreMaster_' kk.tl:0; execute_unload;
* Solve master problem to propose x(kk) consider trust region of current incumbent (k_curinc)
    solve master min zmaster us lp;

    obj_master(kk) = zmaster.l - sum(scen, prob(scen) * theta.l(scen));
    obj_model(kk,scen) = theta.l(scen);

* Save current solution
    p_PLANT_ADD_CAP(kk,rp)    = PLANT_ADD_CAP.l(rp)  ;
    p_STORAGE_ADD_CAP(kk,rs)  = STORAGE_ADD_CAP.l(rs);
    p_LINK_ADD_CAP(kk,net)    = LINK_ADD_CAP.l(net)  ;

    sub_PLANT_ADD_CAP(rp)    = PLANT_ADD_CAP.l(rp)  ;
    sub_STORAGE_ADD_CAP(rs)  = STORAGE_ADD_CAP.l(rs);
    sub_LINK_ADD_CAP(net)    = LINK_ADD_CAP.l(net)  ;


* Fix bounds
*    PLANT_ADD_CAP.fx(rp)    = p_PLANT_ADD_CAP(kk,rp) ;
*    STORAGE_ADD_CAP.fx(rs)  = p_STORAGE_ADD_CAP(kk,rs);
*    LINK_ADD_CAP.fx(net)    = p_LINK_ADD_CAP(kk,net);
    putclose putlog / '##### before Master EC1' /;
    display sub_PLANT_ADD_CAP, sub_STORAGE_ADD_CAP, sub_LINK_ADD_CAP;
    rep(kk,'subTime') = timeelapsed;
    continueEmbeddedCode:
      comm.bcast([[0]] + symbols2List('sub_PLANT_ADD_CAP', 'sub_STORAGE_ADD_CAP', 'sub_LINK_ADD_CAP'), root=0)
      cut = comm.gather(None, root=0)[1:]
      gams.set('sub_obj_sub',      [r for c in cut for r in c[0]])
      gams.set('sub_bcoef_plant' , [r for c in cut for r in c[1]])
      gams.set('sub_bcoef_stor',   [r for c in cut for r in c[2]])
      gams.set('sub_bcoef_link' ,  [r for c in cut for r in c[3]])
      gams.set('sub_bconst',       [r for c in cut for r in c[4]])  
    pauseEmbeddedCode sub_obj_sub sub_bcoef_plant sub_bcoef_stor sub_bcoef_link sub_bconst
	display  sub_obj_sub, sub_bcoef_plant, sub_bcoef_stor, sub_bcoef_link, sub_bconst;
    rep(kk,'subTime') = timeelapsed - rep(kk,'subTime'); 
    putclose putlog / '##### after Master EC1' /;
    
*   copy non-k-indexed sub problem data to k-indexed original data!!!!
    obj_sub(kk,scen)         = sub_obj_sub(scen)       ;
    bcoef_plant(kk,scen,rp)  = sub_bcoef_plant(scen,rp);
    bcoef_stor(kk,scen,rs)   = sub_bcoef_stor(scen,rs) ;
    bcoef_link(kk,scen,net)  = sub_bcoef_link(scen,net);
    bconst(kk,scen)          = sub_bconst(scen)        ; 
    
* After all subproblems have been evaluated activate new cut for next master problem    
    kcl_cut(kk,cl) = yes;

* Reset bounds
    PLANT_ADD_CAP.lo(rp)    = 0;
    STORAGE_ADD_CAP.lo(rs)  = 0;
    LINK_ADD_CAP.lo(net)    = 0;
    PLANT_ADD_CAP.up(rp)    = plant_max_add_cap(rp);
    STORAGE_ADD_CAP.up(rs)  = storage_max_add_cap(rs);
    LINK_ADD_CAP.up(net)    = link_max_add_cap(net);


* if solution performs worse than current incumbent and decrease trust region for large differences between model and sub objective or consecutive errors
    if( (sum(scen, prob(scen) * obj_sub(kk,scen)) + obj_master(kk))  >  (sum( (scen,k1)$k_curinc(k1), prob(scen) * obj_sub(k1,scen) ) + sum( k1$k_curinc(k1), obj_master(k1) )) - eta * ( (sum( (scen,k1)$k_curinc(k1), prob(scen) * obj_sub(k1,scen) ) + sum( k1$k_curinc(k1), obj_master(k1) )) - (sum(scen, prob(scen) * obj_model(kk,scen)) + obj_master(kk)) ) ,

$IFTHEN.useDynamicTrustRegion %DYNAMICTRUSTREG%==1
        trr_rho = min(1,trr(kk)) * ( (sum(scen, prob(scen)*obj_sub(kk,scen))+obj_master(kk)) - (sum((scen,k1)$k_curinc(k1),prob(scen)*obj_sub(k1,scen))+sum(k1$k_curinc(k1),obj_master(k1))) ) / ( (sum((scen,k1)$k_curinc(k1),prob(scen)*obj_sub(k1,scen))+sum(k1$k_curinc(k1),obj_master(k1))) - (sum(scen, prob(scen)*obj_model(kk,scen))+obj_master(kk)) );

        if(trr_rho > 0,
            trr_counter = trr_counter + 1;
        );
        
        if(trr_rho > 3 OR ( trr_counter >= 3 AND (trr_rho > 1 AND trr_rho <= 3) ),

            trr(k)$(ord(k) > ord(kk)) = max(trr_min, 1/min(trr_rho, 4) * trr(k));

            display$(loglevel >= 2) 'Trust region updated',
                    trr;

            trr_counter = 0;
        );
$ENDIF.useDynamicTrustRegion

* ... else increase trust region radius for small differences between model and sub objective...
    else

$IFTHEN.useDynamicTrustRegion %DYNAMICTRUSTREG%==1
        if( (sum(scen, prob(scen) * obj_sub(kk,scen)) + obj_master(kk))  <=  (sum( (scen,k1)$k_curinc(k1), prob(scen) * obj_sub(k1,scen) ) + sum( k1$k_curinc(k1), obj_master(k1) )) - 0.5 * ( (sum( (scen,k1)$k_curinc(k1), prob(scen) * obj_sub(k1,scen) ) + sum( k1$k_curinc(k1), obj_master(k1) )) - (sum(scen, prob(scen) * obj_model(kk,scen)) + obj_master(kk)) ) ,
            trr(k)$(ord(k) > ord(kk)) = min( 2 * trr(kk),trr_max );

            display$(loglevel >= 2) 'Trust region updated',
                    trr;
        );
$ENDIF.useDynamicTrustRegion

* ... then accept current incumbent as new incumbent (necessary for new trust region)
        k_curinc(k) = no;
        k_curinc(kk) = yes;

        display$(loglevel >= 2) 'New incumbent set',
                k_curinc;

    );


* ... check for optimality ...
    upperbound(kk)  = obj_master(kk) + sum(scen, prob(scen) * obj_sub(kk,scen));
    upperbound_glob = min(upperbound_glob, upperbound(kk));
    lowerbound(kk)  = obj_master(kk) + sum(scen, prob(scen) * obj_model(kk,scen));
    
    conv_gap(kk) = upperbound(kk) - lowerbound(kk);
    display$(1) 'Current convergence gap:',
            conv_gap, upperbound, lowerbound;

    obj_total(kk) = obj_master(kk) + sum(scen, prob(scen) * obj_sub(kk,scen));
    display$(0) obj_total, obj_master, obj_model, obj_sub;

* Check convergence
    done = conv_gap(kk) <= %CONVERGENCECRIT%;

    rep(kk,'Obj')     = upperbound(kk);
    rep(kk,'globObj') = upperbound_glob;
    rep(kk,'Bound')   = lowerbound(kk);
    rep(kk,'absGap')  = conv_gap(kk);
    rep(kk,'relGap')  = conv_gap(kk)/upperbound(kk);
*   rep(kk,'subTime') calculated above
    rep(kk,'masterTime') = timeelapsed - rep(kk,'masterTime') - rep(kk,'subTime');    
    rep(kk,'Time')    = timeelapsed;
* Log time of point evaluated
    t_e_k(kk) = timeelapsed;
);
* Terminate sub jobs
continueEmbeddedCode:
  comm.bcast([[1],[],[],[]], root=0)
endEmbeddedCode
* report summary
rep('TOTAL','subTime')    = sum(kk$rep(kk,'subTime'), rep(kk,'subTime'));
rep('TOTAL','masterTime') = sum(kk$rep(kk,'masterTime'), rep(kk,'masterTime'));
rep('TOTAL','Time')       = timeelapsed;

* Get best solution
Parameter
final_PLANT_ADD_CAP(rr,p)
final_STORAGE_ADD_CAP(rr,s)
final_LINK_ADD_CAP(rr1,rr2)
final_obj_total
final_obj_master
final_obj_sub;

final_PLANT_ADD_CAP(rp)   = sum(k$k_curinc(k), p_PLANT_ADD_CAP(k,rp));
final_STORAGE_ADD_CAP(rs) = sum(k$k_curinc(k), p_STORAGE_ADD_CAP(k,rs));
final_LINK_ADD_CAP(net)   = sum(k$k_curinc(k), p_LINK_ADD_CAP(k,net));
final_obj_total           = sum(k$k_curinc(k), obj_total(k));
final_obj_master          = sum(k$k_curinc(k), obj_master(k));
final_obj_sub             = sum(k$k_curinc(k), sum(scen,prob(scen) * obj_sub(k,scen)));

display "-----------------------------------------------------------------",
        "Optimal solution and cost",
        final_PLANT_ADD_CAP,final_STORAGE_ADD_CAP,final_LINK_ADD_CAP,final_obj_total,final_obj_master,final_obj_sub,
        "-----------------------------------------------------------------";
t_e = timeelapsed;
display "-----------------------------------------------------------------",
        " Time elapsed: ", t_e,
        "-----------------------------------------------------------------";
display rep;
Execute_Unload 'spBendersMPI_results_%NBREGIONS%.gdx',done,t_e,t_e_k,trr,clscen,kcl_cut,k_incmap,obj_total,obj_master,obj_model,obj_sub,bconst,bconst_check,bcoef_plant,bcoef_stor,bcoef_link,p_PLANT_ADD_CAP,p_STORAGE_ADD_CAP,p_LINK_ADD_CAP,final_PLANT_ADD_CAP,final_STORAGE_ADD_CAP,final_LINK_ADD_CAP;

$else.MPI
  activeScen(scen) = ord(scen)=%sysenv.PMI_RANK%;

  while(1,
    putclose putlog '##### before Sub EC1';
    continueEmbeddedCode:
    primal_solution = comm.bcast(None, root=0)
    gams.set('done', primal_solution[0])
    gams.set('sub_PLANT_ADD_CAP', primal_solution[1])
    gams.set('sub_STORAGE_ADD_CAP', primal_solution[2])
    gams.set('sub_LINK_ADD_CAP', primal_solution[3])
    pauseEmbeddedCode done sub_PLANT_ADD_CAP sub_STORAGE_ADD_CAP sub_LINK_ADD_CAP 
    abort.noerror$done 'terminating subprocess';
    putclose putlog '##### after Sub EC1';
    display sub_PLANT_ADD_CAP, sub_STORAGE_ADD_CAP, sub_LINK_ADD_CAP;
   
*   Evaluate subgradients at x(kk) for each scenario
    loop(activeScen(scen),

      PLANT_ADD_CAP.fx(rp)  = sub_PLANT_ADD_CAP(rp)  ;
      STORAGE_ADD_CAP.fx(rs)= sub_STORAGE_ADD_CAP(rs);
      LINK_ADD_CAP.fx(net)  = sub_LINK_ADD_CAP(net)  ;
      
*   Get current values for the type_mult(type) parameter for each scenario subproblem
      type_mult(type) = type_mult_scen(scen,type)
      
      solve sub min zsub us lp;
    
      sub_obj_sub(scen) = zsub.l;
    
      sub_bcoef_plant(scen,rp)  = eps + sum(t, eq_plant_capacity.m(t,rp) * avail(t,rp) * %RESOLUTION%);
      sub_bcoef_stor(scen,rs)   = eps + sum(t, eq_storage_capacity.m(t,rs) );
      sub_bcoef_link(scen,net)  = eps + sum(t, eq_link_capacity.m(t,net) * %RESOLUTION% );
    
      sub_bconst(scen) =  eps + sub_obj_sub(scen)
                       - sum(rp,  sub_bcoef_plant(scen,rp) * sub_PLANT_ADD_CAP(rp))
                       - sum(rs,  sub_bcoef_stor(scen,rs)  * sub_STORAGE_ADD_CAP(rs))
                       - sum(net, sub_bcoef_link(scen,net) * sub_LINK_ADD_CAP(net));
    );
    display  sub_bcoef_plant;
    continueEmbeddedCode:
    comm.gather(symbols2List('sub_obj_sub', 'sub_bcoef_plant', 'sub_bcoef_stor', 'sub_bcoef_link', 'sub_bconst'), root=0 )
    pauseEmbeddedCode
    
  );  
$endif.MPI
