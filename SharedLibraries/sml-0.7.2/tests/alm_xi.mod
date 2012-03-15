# This file is (c) 2008 Xi Yang, Andreas Grothey, University of Edinburgh
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty.

#
# ALM model with Stochastic dominance constraints written in SML
#

param T;
set TIME ordered := 0..T;
set NODES;
param Parent{NODES} symbolic; # parent of nodes
param Probs{NODES}; # probability distribution of nodes

set ASSETS;
param Price{ASSETS}; # prices of assets
param Return{ASSETS, NODES}; # returns of assets at each node

set BENCHMARK;
param Vbenchmk{BENCHMARK}; # values of benchmarks
param Bench2nd{BENCHMARK}; # 2nd order SD values of benchmarks
param Gamma; # transaction costs

param InitialWealth;

var slack2{BENCHMARK, TIME} >= 0;

block alm stochastic using (nd in NODES, Parent, Probs, st in TIME): {

  var x_hold{ASSETS} >= 0;
  var risk{BENCHMARK} >= 0;
  var cash >= 0;

  stages {0}: {
    subject to StartBudget:
      (1+Gamma) * sum{j in ASSETS} x_hold[j] * Price[j] + cash = InitialWealth;
  }

  stages (1..T): {

    var x_sold{ASSETS} >= 0;
    var x_bought{ASSETS} >= 0;
    var slack1{BENCHMARK} >= 0;

    subject to CashBalance:
      cash + (1+Gamma) * sum{j in ASSETS} Price[j] * x_bought[j] =
        ancestor(1).cash + (1-Gamma) * sum{j in ASSETS} Price[j] * x_sold[j];

    subject to Inventory{j in ASSETS}:
      x_hold[j] =
        (1+Return[j,nd]) * ancestor(1).x_hold[j] + x_bought[j] - x_sold[j];

    subject to StochasticDominance1{l in BENCHMARK}:
      sum{j in ASSETS}
        Price[j] * x_hold[j] + cash - Vbenchmk[l] = -risk[l] + slack1[l];

    subject to StochasticDominance2{l in BENCHMARK}:
      Exp(risk[l]) = Bench2nd[l] - slack2[l,st];
  }

  stages {T}: {

    var wealth >= 0;

    subject to FinalWealth:
      wealth = sum{j in ASSETS} Price[j] * x_hold[j] + cash;

    maximize objFunc: wealth;
  }
}
