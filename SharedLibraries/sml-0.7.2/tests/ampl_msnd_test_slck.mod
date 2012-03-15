# This file is (c) 2008 Andreas Grothey, University of Edinburgh
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty.

# this is a structured AMPL model of the MSND problem 
#
#   In MSND a network (with installed capacity) is repeated many times
#   every time with one of the ARCS (or NODES) missing. The aim is to
#   install as little as possible of spare capacity on each arc to make
#   sure that the (constant) demand can be routed for every failure
#


#MSND

set NODES ordered;
set ARCS within NODES cross NODES;
set COMM;

param cost{ARCS};
param cap{ARCS};

param src{COMM};
param dest{COMM};
param amount{COMM};

var sparecap{(i1,i2) in ARCS} >=0;

#the basic MCNF block is repeated for every (missing) arc 
block MCNF{(i1,i2) in ARCS}:
  
  set ARCSDIFF := ARCS diff {(i1,i2)};

  # the RouteComm block routes the j'th commodity on the network with arc i
  # missing
  block RouteComm{j in COMM}:
     var Flow{(s,t) in ARCSDIFF} >=0;
     subject to FlowBalance{l in NODES}:
         sum{(m,l) in ARCSDIFF} Flow[m,l] + (if (ord(l)==src[j]) then amount[j]) =
         sum{(l,m) in ARCSDIFF} Flow[l,m] + (if (ord(l)==dest[j]) then amount[j]);
#     subject to FlowBalance{l in NODES}:
#         sum{m in ARCSDIFF} Flow[l] =
#         sum{m in ARCSDIFF} Flow[l,m];
  end block;

  var CapFlowSlack{(k1,k2) in ARCSDIFF}>=0;
  subject to CapFlow{(k1,k2) in ARCSDIFF}:
       sum{j in COMM} (RouteComm[j].Flow[k1,k2]) +CapFlowSlack[k1,k2]= cap[k1,k2] + sparecap[k1,k2];
#       sum{j in COMM} (RouteComm[j].Flow[k].This[i]) <= cap[k] + sparecap[k];

end block;

minimize obj1: sum{(i,j) in ARCS} sparecap[i,j]*cost[i,j];


