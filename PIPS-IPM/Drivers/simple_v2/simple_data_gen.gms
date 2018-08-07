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

Alias (rr,rr1,rr2);

* Subsets and mappings
Set pr(p,rr)                     'plant to region mapping                        '
    rp(rr,p)                     'region to plant mapping                        '
    ptype(rr,p,type)             'plant type mapping                             '
    sr(s,rr)                     'storage to region mapping                      '
    rs(rr,s)                     'region to storage mapping                      '
    net(rr1,rr2)                 'transmission links                             '
    net_l(rr1,rr2)               'transmission links from lower to higher index  '
    net_g(rr1,rr2)               'transmission links from higher to lower index  '
;

Parameter
* time indexed
    plant_capX(rr,p,ttX)         'plant capacity                       [GW]      '
    cost_unserved_demandX(ttX)   'price for unserved electricity load  [MEUR/GWh]'
    link_capX(rr1,rr2,ttX)       'transmission link capacity per hour  [GWh]     '
    link_efficiencyX(rr1,rr2,ttX)'transmission link efficiency factor            '
    demandX(rr,ttX)              'demand for region per time step      [GWh]     '
    ecarsX(rr,ttX)               'demand for ecar energy per region/time step [GWh]     '
    renewablesX(rr,ttX)          'available renewable energy per hour  [GWh]     '
    base_time_series(ttX,*)      'base case time series                [GWh]     '
    base_plant_data(type,*)      'base case plant types                          '
    ecarNorm(ttX)
    demandEcarX(rr,ttX)

* not time indexed
    regionsize(rr)               'region size compared to base case              '
    plantsize(rr,p)              'plant size compared to base case               '
    storagesize(rr,s)            'storage size compared to base case             '
    location(rr,*)               'location of region in a 1000x1000 km square    '
    distance(rr1,rr2)            'distance between regions             [km]      '
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
    link_max_add_cap(rr1,rr2)    'max additional arc capacity per hour [GWh]     '
    cost_plant_add(rr,p)         'cost for additional plant capacity   [MEUR/GW] '
    cost_storage_add(rr,s)       'cost for additional storage capacity [MEUR/GWh]'
    cost_link_add(rr1,rr2)       'cost for additional link capacity    [MEUR/GWh]'
    availibility(ttX,type)       'availability of energy types                   '
    availX(ttX,rr,p)             'availability of plant                          '
    availXrev(rr,p,ttX)          'availability of plant changed index order      '
    availavg(rr,p)               'avarage plant availability                     '
    rpshift(rr,p)                'random shift of time series data for plant     '
    startX(ttX)                  'start hour of data time step                   '
    endX(ttX)                    'end hour of data time step                     '
    demandSum(rr)
;

scalar
    base_plant_cap               'base case plant capacity             [GW]      ' / 75 /
    ecarSum
;

* compute mappings (every region gets a fossil and a renewable plant and one storage)
pr(p,rr) = ord(p)>=2*ord(rr)-1 and ord(p)<=2*ord(rr);
sr(s,rr) = ord(s) = ord(rr);
option rp<pr, rs<sr;

$call csv2gdx %SIMPLEDIR%timeSeries.csv output=%gams.scrdir%timeSeries.gdx ID=base_time_series UseHeader=y Index=1 Values=2..lastCol
$if errorlevel 1 $abort Problems reading %SIMPLEDIR%timeSeries.csv with csv2gdx
$gdxin %gams.scrdir%timeSeries.gdx
$load time_series_head<base_time_series.Dim2 base_time_series
$gdxin

$call csv2gdx %SIMPLEDIR%plantTypeData.csv output=%gams.scrdir%plantTypeData.gdx ID=base_plant_data UseHeader=y Index=1 Values=2..lastCol
$if errorlevel 1 $abort Problems reading %SIMPLEDIR%plantTypeData.csv with csv2gdx
$gdxin %gams.scrdir%plantTypeData.gdx
$load plant_data_head<base_plant_data.Dim2 base_plant_data
$gdxin

* rescale time series and plant data from MWh to GWh and from EUR to MEUR
base_time_series(ttX,time_series_head) = 0.001 * base_time_series(ttX,time_series_head);

* compute availability of renewables
parameter maxre(type);
maxre(re) = smax(ttX, base_time_series(ttX,re));
availibility(ttX,re) = base_time_series(ttX,re)/maxre(re);
availibility(ttX,type)$(not re(type)) = 1;

* random region location in 1000x1000 km roster square
location(rr,'x') = uniform(0,1000);
location(rr,'y') = uniform(0,1000);
distance(rr1,rr2)  = edist(location(rr1,'x')-location(rr2,'x'),location(rr1,'y')-location(rr2,'y'));
* insert links depending on distance between regions an NBNEIGHBOURS
parameter rank(rr1,rr2);

embeddedCode Python:
    rCount = len(gams.get('rr'))
    d = list(gams.get('distance'))
    d = [d[x:x+rCount-1] for x in range(0, len(d), rCount-1)]
    rank = []
    for chunk in d:
        tmp = sorted(chunk, key=lambda x:x[1])
        for idx in range(rCount-1):
            rank.append((tmp[idx][0], idx+1))
    gams.set('rank', rank)
endEmbeddedCode rank
;
net(rr1,rr2) = rank(rr1,rr2) > 0.5 and rank(rr1,rr2) < %NBNEIGHBOURS% + 0.5;

* make net symmetric
net(rr1,rr2) = net(rr1,rr2) or net(rr2,rr1);
net_l(net(rr1,rr2)) = ord(rr1)<ord(rr2);
net_g(net(rr1,rr2)) = ord(rr1)>ord(rr2);

* random plant type
set evenp(p),oddp(p);
evenp(p) = mod(ord(p),2)=0;
oddp(p)  = not evenp(p);
ptype(rp(rr,oddp),'coal')          = uniform(0,1) < 0.5;
ptype(rp(rr,oddp),'gas')           = not ptype(rp,'coal');
ptype(rp(rr,evenp),'photovoltaic') = uniform(0,1) < 0.5;
ptype(rp(rr,evenp),'wind')         = not ptype(rp,'photovoltaic');

availX(ttX,rp(rr,p))               = sum(ptype(rr,p,type), availibility(ttX,type));
* shift wind time series  about 1 - 10 days (24-240 hours)
* shift photovoltaic time series  about 0 - 7 days
rpshift(rp) = uniformint(24,240);
*availX(ttX,rp(rr,p))$ptype(rr,p,'wind')         = availX(ttX++rpshift(rp),rp);
option availXrev<availX;
availXrev(rp(rr,p),ttX)$ptype(rp,'wind')         = availXrev(rp,ttX++rpshift(rp));
rpshift(rp) = uniformint(0,7);
*availX(ttX,rp(rr,p))$ptype(rr,p,'photovoltaic') = availX(ttX++(24*rpshift(rp)),rp);
availXrev(rp(rr,p),ttX)$ptype(rp,'photovoltaic') = availXrev(rp,ttX++(24*rpshift(rp)));
option availX<availXrev;
option clear=availXrev;

availavg(rp) = sum(ttX, availX(ttX,rp)) / card(ttX);

* random plant and region size (base case is germany)
plantsize(rp(rr,p))   = uniform(0.1,1);
regionsize(rr)        = sum(rp(rr,p), plantsize(rp)*availavg(rp)) * uniform(0.9,1.1);
storagesize(rs(rr,s)) = regionsize(rr) * uniform(0.9,1.1);

* populate time indexed parameters
startX(ttX)                        = ord(ttX)-1;
endX(ttX)                          = ord(ttX);
plant_capX(rp,ttX)                 = plantsize(rp) * (base_plant_cap) * uniform(0.95,1.05);
cost_unserved_demandX(ttX)         = 1e4;
link_capX(net_l(rr1,rr2),ttX)      = (sum(rp(rr1,p), plant_capX(rp,ttX)) + sum(rp(rr1,p), plant_capX(rp,ttX))) / 2 * uniform(0.08,0.12);
link_efficiencyX(net_l,ttX)        = max(0.1, 1 - (distance(net_l)/100 * 0.06)) * uniform(0.98,1.02);
demandX(rr,ttX)                    = regionsize(rr) * base_time_series(ttX,'demand') * uniform(0.9,1.0);
ecarSum                            = sum(ttX,base_time_series(ttX,'ecars'));
demandSum(rr)                      = sum(ttX, demandX(rr,ttX))*uniform(0.05,0.15);
ecarNorm(ttX)                      = base_time_series(ttX,'ecars')/ecarSum;
demandEcarX(rr,ttX)                = ecarNorm(ttX)*demandSum(rr)*uniform(0.95,1.05);

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
plant_max_add_cap(rp)              = smax(ttX, plant_capX(rp,ttX)) * uniform(0.2,0.35);
storage_max_add_cap(rs)            = storage_cap(rs) * uniform(0.1,0.25);
link_max_add_cap(net_l)            = sum(ttX,link_capX(net_l,ttX)) / card(ttX) * uniform(0.1,0.5);
cost_plant_add(rp)                 = sum(ptype(rp,type), base_plant_data(type,'expansion_cost')) * uniform(0.9,1.1);;
cost_storage_add(rs)               = 10 * uniform(0.95,1.05);
cost_link_add(net_l)               = uniform(50,200) * 0.4;
* link parameters should be symmetric
link_capX(net_g(rr1,rr2),ttX)       = link_capX(rr2,rr1,ttX) ;
link_efficiencyX(net_g(rr1,rr2),ttX)= link_efficiencyX(rr2,rr1,ttX);
link_max_add_cap(net_g(rr1,rr2))    = link_max_add_cap(rr2,rr1) ;
cost_link_add(net_g(rr1,rr2))       = cost_link_add(rr2,rr1);

$ifthen set WRITEGDX
execute_unload 'simple_data_%NBREGIONS%.gdx';
$endif