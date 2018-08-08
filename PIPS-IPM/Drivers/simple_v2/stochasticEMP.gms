$if not set NBSCEN $set NBSCEN 005
$ifthen.nbscen %NBSCEN%==005
   set scen / scen001*scen005 /;
   parameter
      prob(scen)                 / scen001 0.2, scen002 0.4, scen003 0.2, scen004 0.15, scen005 0.05 /
      scen_type_multX(type,scen) / gas. ( scen001 0.8, scen002 1.0, scen003 1.2, scen004 1.50, scen005 2.00 )
                                   coal.( scen001 1.0, scen002 1.0, scen003 1.1, scen004 1.20, scen005 1.30 )
                                   wind.#scen 1, photovoltaic.#scen 1/
$elseif.nbscen %NBSCEN%==100
   set scen / scen001*scen100 /;
   parameter
      prob(scen)                 / scen001       0.0100  , scen002       0.0100  , scen003       0.0100  , scen004       0.0100  , scen005       0.0100  , scen006       0.0100  , scen007       0.0100  , scen008       0.0100  , scen009       0.0100  , scen010       0.0100  , scen011       0.0100  , scen012       0.0100  , scen013       0.0100  , scen014       0.0100  , scen015       0.0100  , scen016       0.0100  , scen017       0.0100  , scen018       0.0100  , scen019       0.0100  , scen020       0.0100  , scen021       0.0100  , scen022       0.0100  , scen023       0.0100  , scen024       0.0100  , scen025       0.0100  , scen026       0.0100  , scen027       0.0100  , scen028       0.0100  , scen029       0.0100  , scen030       0.0100  , scen031       0.0100  , scen032       0.0100  , scen033       0.0100  , scen034       0.0100  , scen035       0.0100  , scen036       0.0100  , scen037       0.0100  , scen038       0.0100  , scen039       0.0100  , scen040       0.0100  , scen041       0.0100  , scen042       0.0100  , scen043       0.0100  , scen044       0.0100  , scen045       0.0100  , scen046       0.0100  , scen047       0.0100  , scen048       0.0100  , scen049       0.0100  , scen050       0.0100  , scen051       0.0100  , scen052       0.0100  , scen053       0.0100  , scen054       0.0100  , scen055       0.0100  , scen056       0.0100  , scen057       0.0100  , scen058       0.0100  , scen059       0.0100  , scen060       0.0100  , scen061       0.0100  , scen062       0.0100  , scen063       0.0100  , scen064       0.0100  , scen065       0.0100  , scen066       0.0100  , scen067       0.0100  , scen068       0.0100  , scen069       0.0100  , scen070       0.0100  , scen071       0.0100  , scen072       0.0100  , scen073       0.0100  , scen074       0.0100  , scen075       0.0100  , scen076       0.0100  , scen077       0.0100  , scen078       0.0100  , scen079       0.0100  , scen080       0.0100  , scen081       0.0100  , scen082       0.0100  , scen083       0.0100  , scen084       0.0100  , scen085       0.0100  , scen086       0.0100  , scen087       0.0100  , scen088       0.0100  , scen089       0.0100  , scen090       0.0100  , scen091       0.0100  , scen092       0.0100  , scen093       0.0100  , scen094       0.0100  , scen095       0.0100  , scen096       0.0100  , scen097       0.0100  , scen098       0.0100  , scen099       0.0100  , scen100       0.0100  /
      scen_type_multX(type,scen) / gas. ( scen001       0.6140  , scen002       1.8220  , scen003       1.8530  , scen004       0.7200  , scen005       1.0180  , scen006       0.7200  , scen007       1.2050  , scen008       1.6220  , scen009       1.3210  , scen010       2.1890  , scen011       1.5550  , scen012       1.0540  , scen013       2.3530  , scen014       1.5020  , scen015       2.5800  , scen016       1.8000  , scen017       0.8110  , scen018       1.8790  , scen019       2.5750  , scen020       2.2880  , scen021       1.9060  , scen022       2.6080  , scen023       2.3300  , scen024       2.9820  , scen025       0.7620  , scen026       2.8520  , scen027       2.2660  , scen028       1.0120  , scen029       2.1740  , scen030       1.6770  , scen031       1.4550  , scen032       1.0420  , scen033       1.2000  , scen034       2.3000  , scen035       1.4000  , scen036       0.9870  , scen037       2.8050  , scen038       2.9550  , scen039       1.5790  , scen040       0.6270  , scen041       1.5070  , scen042       1.8130  , scen043       2.6840  , scen044       1.0110  , scen045       1.6950  , scen046       0.9740  , scen047       1.2280  , scen048       2.6530  , scen049       1.2010  , scen050       2.5680  , scen051       2.1520  , scen052       0.7310  , scen053       2.2760  , scen054       2.2260  , scen055       2.7570  , scen056       1.8530  , scen057       1.1140  , scen058       1.0240  , scen059       1.9700  , scen060       1.9570  , scen061       1.8150  , scen062       2.2790  , scen063       2.8660  , scen064       1.4630  , scen065       1.1160  , scen066       2.1920  , scen067       2.9350  , scen068       1.4110  , scen069       0.8070  , scen070       2.8570  , scen071       1.3400  , scen072       0.9180  , scen073       1.7590  , scen074       0.7760  , scen075       0.6110  , scen076       2.2480  , scen077       1.7850  , scen078       2.5860  , scen079       1.0770  , scen080       2.7750  , scen081       1.4270  , scen082       2.3250  , scen083       2.5790  , scen084       1.2390  , scen085       2.8460  , scen086       2.5940  , scen087       1.4980  , scen088       1.8790  , scen089       0.7780  , scen090       1.4450  , scen091       1.7530  , scen092       2.7590  , scen093       1.2090  , scen094       2.3950  , scen095       1.4190  , scen096       2.1640  , scen097       0.8880  , scen098       2.6770  , scen099       1.5540  , scen100       1.6770 )
                                   coal.( scen001       0.9910  , scen002       1.0140  , scen003       1.6100  , scen004       2.6430  , scen005       2.0500  , scen006       0.9370  , scen007       1.3000  , scen008       1.1010  , scen009       2.3920  , scen010       0.6100  , scen011       1.5560  , scen012       1.0550  , scen013       2.7040  , scen014       2.6480  , scen015       2.9180  , scen016       1.8820  , scen017       1.0830  , scen018       0.7930  , scen019       0.8580  , scen020       0.9590  , scen021       2.3170  , scen022       1.1280  , scen023       1.0060  , scen024       2.0230  , scen025       1.8500  , scen026       2.7760  , scen027       1.9680  , scen028       0.7380  , scen029       2.2510  , scen030       1.9380  , scen031       2.6860  , scen032       2.1360  , scen033       0.9030  , scen034       1.1620  , scen035       2.1920  , scen036       1.7490  , scen037       2.2830  , scen038       0.7260  , scen039       1.5270  , scen040       2.7960  , scen041       0.7240  , scen042       1.8650  , scen043       1.8880  , scen044       2.3560  , scen045       2.4810  , scen046       2.9040  , scen047       2.9110  , scen048       2.8260  , scen049       2.4360  , scen050       0.7260  , scen051       2.6470  , scen052       1.1660  , scen053       0.9660  , scen054       2.8260  , scen055       1.9950  , scen056       1.8280  , scen057       2.8190  , scen058       2.2120  , scen059       0.8740  , scen060       2.6890  , scen061       1.9630  , scen062       1.0290  , scen063       1.2140  , scen064       2.1610  , scen065       1.3370  , scen066       1.4010  , scen067       2.7080  , scen068       2.7110  , scen069       0.9060  , scen070       1.2430  , scen071       0.9850  , scen072       0.7270  , scen073       1.5800  , scen074       2.1820  , scen075       2.0940  , scen076       2.0120  , scen077       1.3560  , scen078       1.8090  , scen079       2.2950  , scen080       2.4160  , scen081       2.3560  , scen082       2.0130  , scen083       2.3940  , scen084       2.8850  , scen085       2.5630  , scen086       0.6910  , scen087       1.5300  , scen088       1.9130  , scen089       2.0650  , scen090       0.9990  , scen091       2.2190  , scen092       1.1800  , scen093       0.8530  , scen094       2.3620  , scen095       1.3480  , scen096       1.4870  , scen097       1.3950  , scen098       1.5310  , scen099       0.9520  , scen100       1.7490 ) /
$else.nbscen
$  abort invalid number of scenarios
$endif.nbscen
      srep(scen,*)       scenario probability / #scen.prob 0 /
      scen_type_mult(scen,type)
   ;
   option scen_type_mult<scen_type_multX;
   Set dictSP / scen      . scenario . ''
                ''        . opt      . srep
                type_mult . randvar  . scen_type_mult /;

   file emp / '%emp.info%' /; put emp '* problem %gams.i%';

   put / 'jrandvar' ;
   loop(type$(not re(type)), put ' ' type_mult.tn(type));
   loop(scen$prob(scen),
     put '   ' prob(scen);
     loop(type$(not re(type)),
       put ' ' scen_type_mult(scen,type);
     );
   );
   putclose / 'stage 1 LINK_ADD_CAP PLANT_ADD_CAP STORAGE_ADD_CAP                                                                '
            / 'stage 2 FLOW POWER SLACK STORAGE_LEVEL STORAGE_INFLOW STORAGE_OUTFLOW EMISSION_SPLIT EMISSION_COST ROBJ OBJ       '
            / 'stage 2 eq_robj eq_power_balance eq_plant_capacity eq_total_plant_capacity eq_storage_balance eq_storage_capacity '
            / 'stage 2 eq_emission_region eq_emission_cost eq_emission_cap eq_link_capacity eq_obj                               '
            / 'stage 2 type_mult'
   ;
   t(tt) = yes;
   r(rr) = yes;

$IFTHENI %EMPMETHOD%==DE
   option emp = de;
$ELSEIFI %EMPMETHOD%==LINDO
   option emp = lindo;
$ELSEIFI %EMPMETHOD%==LINDOBENDERS
   option emp = lindo;
   simple.optfile = 1;
$  onecho  > lindo.opt
STOC_DEQOPT     10
STOC_MAP_MPI2LP  1
$  offecho
$ELSE
$  abort 'invalid setting of EMPMETHOD. Allowed are "DE", "LINDO", "LINDOBENDERS".
$ENDIF
   solve simple min OBJ use emp scenario dictSP;
