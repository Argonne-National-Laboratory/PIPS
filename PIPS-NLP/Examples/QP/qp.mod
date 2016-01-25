param NS;
param NS_ALL;
set SCEN;                                # scenario set
set SCEN_ALL:=1..NS_ALL;                 # set of all scenarios

param nz;
set  FirstVar:=1..nz;                
param zl{FirstVar}; 
param zu{FirstVar}; 
param dz{FirstVar};
param Hz{FirstVar,FirstVar};

param nw;
set  SecondVar:=1..nw;
param wl{SecondVar};
param wu{SecondVar};
param dw{SecondVar};
param Hw{SecondVar,SecondVar};

param nec;
set EConstraint:=1..nec;
param bEw{SCEN_ALL, EConstraint};
param TEz{EConstraint,FirstVar};
param TEw{EConstraint,SecondVar};

param nic;
set IConstraint:=1..nic;
param cl{SCEN_ALL, IConstraint};
param cu{SCEN_ALL, IConstraint};
param TIz{IConstraint,FirstVar};
param TIw{IConstraint,SecondVar};


var z{ k in FirstVar}   >= zl[k], <=zu[k], := 0.5;
var w{ s in SCEN, k in SecondVar}   >= wl[k], <=wu[k], := 0.6;


minimize obj: (sum{k in FirstVar} z[k]*dz[k] + 0.5* sum{k in FirstVar, j in FirstVar} z[k]*Hz[k,j]*z[j]   + 1.0/NS*sum{s in SCEN, k in SecondVar} w[s,k]*dw[k]+ 0.5/NS * sum{s in SCEN, k in SecondVar, j in SecondVar} w[s,k]*Hw[k,j]*w[s,j]) ;


Econ{s in SCEN, i in EConstraint}: sum{k in FirstVar}TEz[i,k]* z[k] +  sum{k in SecondVar}TEw[i,k]* w[s,k]           = bEw[s,i];

Icon{s in SCEN, i in IConstraint}: cl[s,i]<= sum{k in FirstVar}TIz[i,k]* z[k] +  sum{k in SecondVar}TIw[i,k]* w[s,k] <= cu[s,i];
