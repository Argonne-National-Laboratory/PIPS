# This file is (c) 2009 Andreas Grothey, Marco Colombo, University of Edinburgh.
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty.

#
# Basic ALM model written in SML
#

param T;
set TIME ordered := 0..T;
set NODES;
param Parent{NODES} symbolic;
param Probs{NODES};

set ASSETS;
param Price{ASSETS};
param Return{ASSETS, NODES};
param Liability{TIME};

param InitialWealth;
param Gamma;  # transaction costs

var makemodelrun >= 0; # at the moment we must have at least one variable here

block alm stochastic using (nd in NODES, Parent, Probs, st in TIME): {

  var x_hold{ASSETS} >= 0;
  var cash >= 0;

  stages {0}: {
    subject to StartBudget:
      (1+Gamma) * sum{j in ASSETS} x_hold[j] * Price[j] + cash = InitialWealth;
  }

  stages (1..T): {

    var x_bought{ASSETS} >= 0;
    var x_sold{ASSETS} >= 0;

    subject to Inventory{j in ASSETS}:
      x_hold[j] =
        (1+Return[j,nd]) * ancestor(1).x_hold[j] + x_bought[j] - x_sold[j];

    subject to CashBalance:
      Liability[st]+cash + (1+Gamma) * sum{j in ASSETS} Price[j] * x_bought[j] =
        ancestor(1).cash + (1-Gamma) * sum{j in ASSETS} Price[j] * x_sold[j];
  }

  stages {T}: {

    var wealth >= 0;

    subject to FinalWealth:
      wealth = sum{j in ASSETS} Price[j] * x_hold[j] + cash;

    maximize objFunc: wealth;
  }
}
