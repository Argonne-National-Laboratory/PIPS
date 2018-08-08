$if not set CPLEXBENDERS $set CPLEXBENDERS 1
$if not set SCENBLOCK    $set SCENBLOCK   -2


$ifthen.nbscen %NBSCEN%==5
   parameter
      prob(scen)                 / scen1 0.2, scen2 0.4, scen3 0.2, scen4 0.15, scen5 0.05 /
      scen_type_multX(type,scen) / gas. ( scen1 0.8, scen2 1.0, scen3 1.2, scen4 1.50, scen5 2.00 )
                                   coal.( scen1 1.0, scen2 1.0, scen3 1.1, scen4 1.20, scen5 1.30 ) /
   option type_mult<scen_type_multX;
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
   option type_mult<scen_type_multX;
   prob(scen) = 1/card(scen);
$else.nbscen
   parameter
      prob(scen)
      type_mult(scen,type);

* Do it for all scenarios, since we want the same random numbers in every block
   type_mult(scen,type)$(not re(type)) = uniform(0.5,2.5);
   prob(scen) = 1/card(scen);

$endif.nbscen

   t(tt) = yes;
   r(rr) = yes;

$ifthen.cpxbenders not %CPLEXBENDERS%==0
   file fcpx / cplexd.opt /;
   put fcpx 'BendersStrategy 1'
          / 'variables.BendersPartition 0';
   loop(scen,
      put / "FLOW.BendersPartition('" scen.tl:0 "',*,*,*) " ord(scen):0:0
          / "POWER.BendersPartition('" scen.tl:0 "',*,*,*) " ord(scen):0:0
          / "SLACK.BendersPartition('" scen.tl:0 "',*,*) " ord(scen):0:0
          / "STORAGE_LEVEL.BendersPartition('" scen.tl:0 "',*,*,*) " ord(scen):0:0
          / "STORAGE_INFLOW.BendersPartition('" scen.tl:0 "',*,*,*) " ord(scen):0:0
          / "STORAGE_OUTFLOW.BendersPartition('" scen.tl:0 "',*,*,*) " ord(scen):0:0
          / "EMISSION_SPLIT.BendersPartition('" scen.tl:0 "',*,*) " ord(scen):0:0
          / "EMISSION_COST.BendersPartition('" scen.tl:0 "',*,*) " ord(scen):0:0
          / "ROBJ.BendersPartition('" scen.tl:0 "',*) " ord(scen):0:0
   );
   putclose;
   sc(scen) = yes;
   STORAGE_INFLOW.up(sc,tt,rs)  = storage_max_in(rs)     ;
   STORAGE_OUTFLOW.up(sc,tt,rs) = storage_max_out(rs)    ;
   option lp=cplexd; simple.optfile=%CPLEXBENDERS%;
   solve simple min OBJ use lp;
$onecho > cplexd.999
lpmethod 4
solutiontype 2
$offecho
$else.cpxbenders
* First stage variables
   LINK_ADD_CAP.stage(net) = 1;
   PLANT_ADD_CAP.stage(rp(rr,p)) = 1;
   STORAGE_ADD_CAP.stage(rs(rr,s)) = 1;
   eq_obj.stage = card(scen)+2;

   file fopt / "%gams.optdir%convertd.opt" /;
   option lp=convertd; simple.optfile=1;

   if (%SCENBLOCK%=0 or %SCENBLOCK%=-1,
     putclose fopt 'jacobian block0.gdx';
     option clear=sc;
     solve simple min OBJ use lp;
   );

* Second stage variables in block scen
$ifthen set SCENLISTSTART
   set sList(scen) / scen%SCENLISTSTART%*scen%SCENLISTEND% /;
   loop(sList(scen),
$else
   loop(scen$(ord(scen)=%SCENBLOCK% or %SCENBLOCK%<0),
$endif
      FLOW.stage(scen,t,net) = ord(scen)+1;
      POWER.stage(scen,t,rp(rr,p)) = ord(scen)+1;
      SLACK.stage(scen,t,rr) = ord(scen)+1;
      STORAGE_LEVEL.stage(scen,t,rs(rr,s)) = ord(scen)+1;
      STORAGE_INFLOW.stage(scen,t,rs(rr,s)) = ord(scen)+1;
      STORAGE_OUTFLOW.stage(scen,t,rs(rr,s)) = ord(scen)+1;
      EMISSION_SPLIT.stage(scen,rr,e) = ord(scen)+1;
      EMISSION_COST.stage(scen,rr,e) = ord(scen)+1;
      ROBJ.stage(scen,rr) = ord(scen)+1;
* Second stage equations in block scen
      eq_robj.stage(scen,rr) = ord(scen)+1;
      eq_power_balance.stage(scen,t,rr) = ord(scen)+1;
      eq_plant_capacity.stage(scen,t,rp(rr,p)) = ord(scen)+1;
      eq_total_plant_capacity.stage(scen,rp(rr,p)) = ord(scen)+1;
      eq_storage_balance.stage(scen,t,rs(rr,s)) = ord(scen)+1;
      eq_storage_capacity.stage(scen,t,rs(rr,s)) = ord(scen)+1;
      eq_emission_region.stage(scen,rr,e) = ord(scen)+1;
      eq_emission_cost.stage(scen,rr,e) = ord(scen)+1;
      eq_emission_cap.stage(scen,e) = ord(scen)+1;
      eq_link_capacity.stage(scen,t,net) = ord(scen)+1;
* Solve
      option clear=sc; sc(scen)=yes;
      STORAGE_INFLOW.up(sc,tt,rs)  = storage_max_in(rs)     ;
      STORAGE_OUTFLOW.up(sc,tt,rs) = storage_max_out(rs)    ;
      if (%SCENBLOCK%>-2,
        putclose fopt 'jacobian block' ord(scen):0:0 '.gdx';
        solve simple min OBJ use lp;
      );
   );
   if (%SCENBLOCK%=-2,
     putclose fopt 'jacobian blockall.gdx';
     sc(scen) = yes;
     solve simple min OBJ use lp;
   );

$endif.cpxbenders

