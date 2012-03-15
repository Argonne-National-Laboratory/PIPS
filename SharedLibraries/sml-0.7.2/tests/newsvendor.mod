# This file is (c) 2009 Marco Colombo, University of Edinburgh.
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty.

#
# Newsvendor model written in SML
#

set TIME ordered := 1..2;
set NODES;
param Parent{NODES} symbolic;
param Probs{NODES};

param SellingPrice;
param InitBuyCost;
param ExtraBuyCost;
param Demand{NODES};
param MaxInitBuy;

var initBuy >= 0;

block sml stochastic using (nd in NODES, Parent, Probs, TIME): {

  stages {1}: {

    var slack >= 0;

    subject to InitialBuyingLimit:
      initBuy + slack = MaxInitBuy;
  }

  stages {2}: {

    var extraBuy >= 0;
    var unsold >= 0;
    var profit;

    subject to DemandSatisfaction:
      initBuy + extraBuy = Demand[nd] + unsold;

    subject to ProfitComputation:
      profit = Demand[nd] * SellingPrice - initBuy * InitBuyCost - extraBuy * ExtraBuyCost;

    maximize Objective:
      profit;
  }
}
