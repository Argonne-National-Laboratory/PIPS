* make up data for one year

$if not set NBREGIONS    $set  NBREGIONS    4
$if not set NBNEIGHBOURS $eval NBNEIGHBOURS min(floor(%NBREGIONS%/2),4)
$eval NBPLANTS %NBREGIONS%*2
$setnames "%gams.input%" SIMPLEDIR fn fe

* Basic sets
Set rr                       regions                    / r1*r%NBREGIONS% /
    p                        plants                     / p1*p%NBPLANTS% /
    s                        storages                   / s1*s%NBREGIONS% /
    ttX                      time steps                 / tX1*tX8760 /
    e                        emissions                  / CO2 /
    type                     plant type                 / coal, gas, photovoltaic, wind /
    re(type)                 renewable energy types     / photovoltaic, wind /
    time_series_head(*)      time series header
    plant_data_head(*)       plant data header
;

Alias (rr,r1,r2);

* Subsets and mappings
Set pr(p,rr)                     'plant to region mapping                        '
    rp(rr,p)                     'region to plant mapping                        '
    ptype(rr,p,type)             'plant type mapping                             '
    sr(s,rr)                     'storage to region mapping                      '
    rs(rr,s)                     'region to storage mapping                      '
    net(r1,r2)                   'transmission links                             '
    net_l(r1,r2)                 'transmission links from lower to higher index  '
    net_g(r1,r2)                 'transmission links from higher to lower index  '
;

Parameter
* time indexed
    plant_capX(rr,p,ttX)         'plant capacity                       [GW]      '
    cost_unserved_demandX(ttX)   'price for unserved electricity load  [MEUR/GWh]'
    link_capX(r1,r2,ttX)         'transmission link capacity per hour  [GWh]     '
    link_efficiencyX(r1,r2,ttX)  'transmission link efficiency factor            '
    demandX(rr,ttX)              'demand for region per time step      [GWh]     '
    renewablesX(rr,ttX)          'available renewable energy per hour  [GWh]     '
    base_time_series(ttX,*)      'base case time series                [GWh]     '
    base_plant_data(type,*)      'base case plant types                          '
* not time indexed
    regionsize(rr)               'region size compared to base case              '
    plantsize(rr,p)              'plant size compared to base case               '
    storagesize(rr,s)            'storage size compared to base case             '
    location(rr,*)               'location of region in a 1000x1000 km square    '
    distance(r1,r2)              'distance between regions             [km]      '
    yearly_plant_cap(rr,p)       'yearly plant capacity                [GWh]     '
    plant_emission(rr,p,e)       'plant emission                       [tons/GWh]'
    cost_power_generation(rr,p)  'electricity production cost          [MEUR/GWh]'
    yearly_emission_cap(e)       'emission cap for time horizon        [tons]    '
    cost_emission(e)             'emission costs                       [MEUR/ton]'
    storage_cap(rr,s)            'storage capacity                     [GWh]     '
    storage_efficiency(rr,s)     'effieciency factor of storage                  '
    storage_efficiency_in(rr,s)  'effieciency factor of storage inflow           '
    storage_efficiency_out(rr,s) 'effieciency factor of storage outflow          '
    storage_max_in(rr,s)         'maximum inflow into storage per hour [GWh]     '
    storage_max_out(rr,s)        'maximum outflow into storage per hour[GWh]     '
    plant_max_add_cap(rr,p)      'max additional plant capacity        [GW]      '
    storage_max_add_cap(rr,s)    'max additional storage capacity      [GWh]     '
    link_max_add_cap(r1,r2)      'max additional arc capacity per hour [GWh]     '
    cost_plant_add(rr,p)         'cost for additional plant capacity   [MEUR/GW] '
    cost_storage_add(rr,s)       'cost for additional storage capacity [MEUR/GWh]'
    cost_link_add(r1,r2)         'cost for additional link capacity    [MEUR/GWh]'
    avail(ttX,type)              'availability of energy types                   '
    availX(ttX,rr,p)             'availability of plant                          '
    availavg(rr,p)               'avarage plant availability                     '
;

scalar base_plant_cap            'base case plant capacity             [GW]      ' / 70 /
;

* compute mappings (every region gets a fossil and a renewable plant and one storage)
pr(p,rr) = ord(p)>=2*ord(rr)-1 and ord(p)<=2*ord(rr);
sr(s,rr) = ord(s) = ord(rr);
option rp<pr, rs<sr;

* read time base case data
$ifthen not exist "%SIMPLEDIR%data.gdx"
$  ifi %system.fileSys%==UNIX $abort No data.gdx and no way to run gdxxrw.
$  call msappavail -Excel
$  if errorlevel 1 $abort No data.gdx and Excel not available.
$  onecho > instructions.txt
dset = time_series_head  rng=timeSeries!A1    rdim=0 cdim=1
par  = base_time_series  rng=timeSeries!A1    rdim=1 cdim=1
dset = plant_data_head   rng=plantTypeData!A1 rdim=0 cdim=1
par  = base_plant_data   rng=plantTypeData!A1 rdim=1 cdim=1
$  offecho
$  call gdxxrw "%SIMPLEDIR%data.xlsx" o="%SIMPLEDIR%data.gdx" @instructions
$  if errorlevel 1 $abort Problems running gdxxrw.
$endif

$gdxin "%SIMPLEDIR%data.gdx"
$load time_series_head base_time_series plant_data_head base_plant_data
$gdxin

* rescale time series and plant data from MWh to GWh and from EUR to MEUR
base_time_series(ttX,time_series_head) = 0.001 * base_time_series(ttX,time_series_head);

* compute availability of renewables
parameter maxre(type);
maxre(re) = smax(ttX, base_time_series(ttX,re));
avail(ttX,re) = base_time_series(ttX,re)/maxre(re);
avail(ttX,type)$(not re(type)) = 1;

* random region location in 1000x1000 km roster square
location(rr,'x') = uniform(0,1000);
location(rr,'y') = uniform(0,1000);
distance(r1,r2)  = edist(location(r1,'x')-location(r2,'x'),location(r1,'y')-location(r2,'y'));
* insert links depending on distance between regions an NBNEIGHBOURS
parameter unsorted(rr), rank(rr);
$onecho > sort.gms
set       rr;
parameter unsorted(rr)
          rank(rr);
$gdxin sortin.gdx
$load rr unsorted
$gdxin
$ LIBINCLUDE rank unsorted rr rank
execute_unload 'sortout.gdx', rank;
$offecho
loop(r1,
  unsorted(r2) = distance(r1,r2);
  execute_unload 'sortin.gdx', rr unsorted;
  execute 'gams sort.gms lo=2';
  execute_load 'sortout.gdx', rank;
  net(r1,r2) = rank(r2) > 1.5 and rank(r2) < %NBNEIGHBOURS% + 1.5;
);
* make net symmetric
net(r1,r2) = net(r1,r2) or net(r2,r1);
net_l(net(r1,r2)) = ord(r1)<ord(r2);
net_g(net(r1,r2)) = ord(r1)>ord(r2);

* random plant type
set evenp(p),oddp(p);
evenp(p) = mod(ord(p),2)=0;
oddp(p)  = not evenp(p);
ptype(rp(rr,oddp),'coal')          = uniform(0,1) < 0.5;
ptype(rp(rr,oddp),'gas')           = not ptype(rp,'coal');
ptype(rp(rr,evenp),'photovoltaic') = uniform(0,1) < 0.5;
ptype(rp(rr,evenp),'wind')         = not ptype(rp,'photovoltaic');

availX(ttX,rp(rr,p))               = sum(ptype(rr,p,type), avail(ttX,type));
* shift wind time series  about 100 - 200 hours
* shift photovoltaic time series  about 0 - 7 days
availX(ttX,rp(rr,p))$ptype(rr,p,'wind')         = availX(ttX++uniformint(24,240),rp);
availX(ttX,rp(rr,p))$ptype(rr,p,'photovoltaic') = availX(ttX++(24*uniformint(0,7)),rp);
availavg(rp) = sum(ttX, availX(ttX,rp)) / card(ttX);

* random plant and region size (base case is germany)
plantsize(rp(rr,p))   = uniform(0.1,1);
regionsize(rr)        = sum(rp(rr,p), plantsize(rp)*availavg(rp)) * uniform(0.9,1.1);
storagesize(rs(rr,s)) = regionsize(rr) * uniform(0.9,1.1);

* populate time indexed parameters
plant_capX(rp,ttX)                 = plantsize(rp) * (base_plant_cap) * uniform(0.95,1.05);
cost_unserved_demandX(ttX)         = 1e4;
link_capX(net_l(r1,r2),ttX)        = (sum(rp(r1,p), plant_capX(rp,ttX)) + sum(rp(r1,p), plant_capX(rp,ttX))) / 2 * uniform(0.08,0.12);
link_efficiencyX(net_l,ttX)        = max(0.1, 1 - (distance(net_l)/100 * 0.06)) * uniform(0.98,1.02);
demandX(rr,ttX)                    = regionsize(rr) * base_time_series(ttX,'demand') * uniform(0.9,1.0);
*renewablesX(rr,ttX)                = regionsize(rr) * (base_time_series(ttX,'wind') + base_time_series(ttX,'photovoltaic')) * uniform(0.9,1.1);
* populate parameters without time index
yearly_plant_cap(rp)               = sum(ttX, plant_capX(rp,ttX)) * uniform(0.9,0.96);
plant_emission(rp,'CO2')           = sum(ptype(rp,type), base_plant_data(type,'CO2')) * uniform(0.9,1.1);
cost_power_generation(rp)          = sum(ptype(rp,type), base_plant_data(type,'cost')) * uniform(0.9,1.1);
yearly_emission_cap(e)             = uniform(0.8,0.9) * sum((rp,ttX), plant_capX(rp,ttX)*plant_emission(rp,e));
cost_emission('CO2')               = 10*1e-6;
storage_max_in(rs)                 = storagesize(rs) * 0.5 * uniform(0.9,1.1);
storage_max_out(rs)                = storagesize(rs) * 0.5 * uniform(0.9,1.1);
storage_cap(rs)                    = 20 * storage_max_in(rs) * uniform(0.75,1.25);
storage_efficiency(rs)             = uniform(0.99,1);
storage_efficiency_in(rs)          = uniform(0.88,0.9);
storage_efficiency_out(rs)         = uniform(0.89,0.91);
plant_max_add_cap(rp)              = smax(ttX, plant_capX(rp,ttX)) * uniform(0.1,0.25);
storage_max_add_cap(rs)            = storage_cap(rs) * uniform(0.1,0.25);
link_max_add_cap(net_l)            = sum(ttX,link_capX(net_l,ttX)) / card(ttX) * uniform(0.1,0.5);
cost_plant_add(rp)                 = sum(ptype(rp,type), base_plant_data(type,'expansion_cost')) * uniform(0.9,1.1);;
cost_storage_add(rs)               = 10 * uniform(0.95,1.05);
cost_link_add(net_l)               = uniform(50,200) * 0.4;
* link parameters should be symmetric
link_capX(net_g(r1,r2),ttX)       = link_capX(r2,r1,ttX) ;
link_efficiencyX(net_g(r1,r2),ttX)= link_efficiencyX(r2,r1,ttX);
link_max_add_cap(net_g(r1,r2))    = link_max_add_cap(r2,r1) ;
cost_link_add(net_g(r1,r2))       = cost_link_add(r2,r1);

execute_unload '%SIMPLEDIR%simple_data_%NBREGIONS%.gdx';
